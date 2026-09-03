#!/usr/bin/env python3
"""Copy top GWAS rows (by qvalue_bh) into paths expected by notebooks."""

import csv
import shutil
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
FDR_DIR = PROJECT_ROOT / "outputs" / "gwas" / "post_gwas_fdr"
TOP_N = 150

# (source under post_gwas_fdr, dest dir under outputs/gwas, outfile name)
EXPORTS = [
    ("gene_binary__pyseer_binary_fdr.tsv", "pyseer_gene", "top_hits.tsv"),
    ("defense_binary__pyseer_defense_binary_fdr.tsv", "pyseer_defense", "top_hits.tsv"),
    ("prophage_binary__pyseer_prophage_binary_fdr.tsv", "pyseer_prophage", "top_hits.tsv"),
    ("hostrec_binary__pyseer_host_receptors_binary_fdr.tsv", "pyseer_host_receptors", "top_hits.tsv"),
]


def main() -> int:
    if not FDR_DIR.is_dir():
        print(f"Missing {FDR_DIR} — run on HPC: bash scripts/run_post_gwas_fdr.sh", file=__import__("sys").stderr)
        return 1
    for src_name, subdir, out_name in EXPORTS:
        src = FDR_DIR / src_name
        dest_dir = PROJECT_ROOT / "outputs" / "gwas" / subdir
        dest_dir.mkdir(parents=True, exist_ok=True)
        dest = dest_dir / out_name
        if not src.is_file():
            print(f"skip missing: {src}")
            continue
        rows = []
        with src.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = reader.fieldnames
            if not fields:
                continue
            for row in reader:
                rows.append(row)
        rows = rows[:TOP_N]
        with dest.open("w", newline="") as handle:
            w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        print(f"Wrote {dest} ({len(rows)} rows from {src.name})")
    # Optional: stash full binary FDR for gene as secondary reference
    gene_full = FDR_DIR / "gene_binary__pyseer_binary_fdr.tsv"
    if gene_full.is_file():
        shutil.copy2(gene_full, PROJECT_ROOT / "outputs" / "gwas" / "pyseer_gene" / "top_hits_binary_full_fdr.tsv")
        print("Copied pyseer_gene/top_hits_binary_full_fdr.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
