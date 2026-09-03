#!/usr/bin/env python3
"""
Read one or more gdtools ANNOTATE TSV files, add mutant_id and normalize columns,
then write a plot-ready table (event-level or gene-level collapse).
"""
import argparse
import csv
import sys
from pathlib import Path


def read_tsv(path, mutant_id=None):
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        if not rows:
            return []
        if mutant_id is None:
            mutant_id = Path(path).stem
        for r in rows:
            r["mutant_id"] = mutant_id
        return rows


def main():
    p = argparse.ArgumentParser(description="Normalize mutation tables for Figure 5a")
    p.add_argument("tsvs", nargs="+", help="Per-mutant gdtools TSV files")
    p.add_argument("-o", "--output", required=True, help="Output plot-ready TSV")
    p.add_argument("--mode", choices=["event", "gene"], default="event",
                   help="event = one row per mutation; gene = collapse by gene")
    args = p.parse_args()

    all_rows = []
    for tpath in args.tsvs:
        mid = Path(tpath).stem
        rows = read_tsv(tpath, mutant_id=mid)
        all_rows.extend(rows)

    if not all_rows:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as f:
            f.write("mutant_id\tgene\tmutation_type\tproduct\tnucleotide_change\tamino_acid_change\tpredicted_effect\tkeep_for_plot\n")
        print("No mutations; wrote empty table.", file=sys.stderr)
        return

    # Normalize column names (gdtools may use different headers)
    key_map = {
        "Gene": "gene",
        "gene": "gene",
        "Mutation": "mutation_type",
        "Mutation type": "mutation_type",
        "Type": "mutation_type",
        "Product": "product",
        "product": "product",
        "Nucleotide change": "nucleotide_change",
        "Amino acid change": "amino_acid_change",
        "Predicted effect": "predicted_effect",
        "Position": "position",
        "Seq ID": "seq_id",
    }
    out_cols = ["mutant_id", "gene", "mutation_type", "product", "nucleotide_change",
                "amino_acid_change", "predicted_effect", "seq_id", "position", "keep_for_plot"]
    norm = []
    for r in all_rows:
        row = {"mutant_id": r.get("mutant_id", "")}
        for k, v in r.items():
            kclean = k.strip() if k else ""
            if kclean in key_map:
                row[key_map[kclean]] = v
        for c in out_cols:
            if c not in row:
                row[c] = ""
        row["keep_for_plot"] = "1"  # default keep
        norm.append(row)

    if args.mode == "gene":
        # Collapse by mutant_id + gene: concatenate mutation types / AA changes
        by_key = {}
        for r in norm:
            g = (r.get("gene") or "").strip() or "unknown"
            key = (r["mutant_id"], g)
            if key not in by_key:
                by_key[key] = {
                    "mutant_id": r["mutant_id"],
                    "gene": g,
                    "mutation_type": set(),
                    "product": r.get("product") or "",
                    "nucleotide_change": [],
                    "amino_acid_change": [],
                    "predicted_effect": set(),
                }
            by_key[key]["mutation_type"].add((r.get("mutation_type") or "").strip())
            if r.get("nucleotide_change"):
                by_key[key]["nucleotide_change"].append(r["nucleotide_change"])
            if r.get("amino_acid_change"):
                by_key[key]["amino_acid_change"].append(r["amino_acid_change"])
            if r.get("predicted_effect"):
                by_key[key]["predicted_effect"].add(r["predicted_effect"])
        norm = []
        for (_, _), v in sorted(by_key.items()):
            norm.append({
                "mutant_id": v["mutant_id"],
                "gene": v["gene"],
                "mutation_type": ";".join(sorted(v["mutation_type"]) or ["?"]),
                "product": v["product"],
                "nucleotide_change": ";".join(v["nucleotide_change"][:5]),
                "amino_acid_change": ";".join(v["amino_acid_change"][:5]),
                "predicted_effect": ";".join(sorted(v["predicted_effect"]) or ["?"]),
                "seq_id": "",
                "position": "",
                "keep_for_plot": "1",
            })

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(norm)
    print(f"Wrote {len(norm)} rows to {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
