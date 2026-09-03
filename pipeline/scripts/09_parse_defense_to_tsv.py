#!/usr/bin/env python3
"""Parse per-sample DefenseFinder outputs into defense_genes.tsv.

Reads <sample_dir>/*_defense_finder_genes.tsv (e.g. BOR-01_defense_finder_genes.tsv
from protein input, or contigs_defense_finder_genes.tsv from nucleotide-only runs).
"""

import csv
from pathlib import Path
from typing import Dict, List, Optional

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFENSE_RAW = PROJECT_ROOT / "outputs" / "defense" / "raw"
OUT_TSV = PROJECT_ROOT / "outputs" / "defense" / "defense_genes.tsv"

HEADER = [
    "sample",
    "contig",
    "start",
    "end",
    "strand",
    "locus_tag",
    "defense_family",
    "defense_system",
]


def find_genes_tsv(sample_dir: Path, sample: str) -> Optional[Path]:
    preferred = sample_dir / f"{sample}_defense_finder_genes.tsv"
    if preferred.is_file():
        return preferred
    alt = sample_dir / "contigs_defense_finder_genes.tsv"
    if alt.is_file():
        return alt
    matches = sorted(sample_dir.glob("*_defense_finder_genes.tsv"))
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        # Prefer sample-prefixed if ambiguous
        for m in matches:
            if m.name.startswith(sample + "_"):
                return m
        return matches[0]
    return None


def row_from_defense_finder(sample: str, r: Dict[str, str]) -> Dict[str, str]:
    """Map DefenseFinder gene TSV columns to pipeline schema."""
    contig = (r.get("replicon") or r.get("contig") or "").strip()
    start = (r.get("hit_begin_match") or r.get("start") or "").strip()
    end = (r.get("hit_end_match") or r.get("end") or "").strip()
    strand = (r.get("strand") or "").strip()
    locus_tag = (r.get("hit_id") or r.get("locus_tag") or "").strip()
    defense_family = (r.get("gene_name") or r.get("defense_family") or "").strip()
    defense_system = (
        r.get("sys_id") or r.get("defense_system") or r.get("model_fqn") or ""
    ).strip()
    return {
        "sample": sample,
        "contig": contig,
        "start": start,
        "end": end,
        "strand": strand,
        "locus_tag": locus_tag,
        "defense_family": defense_family,
        "defense_system": defense_system,
    }


def parse_sample(sample_dir: Path) -> List[Dict[str, str]]:
    sample = sample_dir.name
    genes_path = find_genes_tsv(sample_dir, sample)
    if genes_path is None:
        return []
    rows: List[Dict[str, str]] = []
    with genes_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for r in reader:
            if not r:
                continue
            # Skip comment-only pseudo-rows
            if not any((v or "").strip() for v in r.values()):
                continue
            rows.append(row_from_defense_finder(sample, r))
    return rows


def main() -> None:
    all_rows: List[Dict[str, str]] = []
    if not DEFENSE_RAW.is_dir():
        raise FileNotFoundError(f"Missing {DEFENSE_RAW}")
    for sample_dir in sorted(DEFENSE_RAW.glob("*")):
        if sample_dir.is_dir():
            all_rows.extend(parse_sample(sample_dir))

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_TSV.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"Wrote {OUT_TSV} ({len(all_rows)} rows)")


if __name__ == "__main__":
    main()
