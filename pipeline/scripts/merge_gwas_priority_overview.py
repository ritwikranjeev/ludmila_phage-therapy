#!/usr/bin/env python3
"""Merge top rows from each post_gwas_fdr/*.tsv into one overview table."""

import argparse
import csv
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--indir",
        type=Path,
        default=Path("outputs/gwas/post_gwas_fdr"),
    )
    ap.add_argument(
        "--outfile",
        type=Path,
        default=Path("outputs/gwas/gwas_priority_overview.tsv"),
    )
    ap.add_argument("--per-file", type=int, default=15)
    args = ap.parse_args()

    if not args.indir.is_dir():
        print(f"missing {args.indir}", file=__import__("sys").stderr)
        return 1

    paths = sorted(args.indir.glob("*_fdr.tsv"))
    if not paths:
        print(f"no *_fdr.tsv in {args.indir}", file=__import__("sys").stderr)
        return 1

    all_rows = []
    fieldnames = None
    for path in paths:
        tag = path.name.replace("_fdr.tsv", "")
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                continue
            base_fields = ["gwas_source"] + list(reader.fieldnames)
            if fieldnames is None:
                fieldnames = base_fields
            for i, row in enumerate(reader):
                if i >= args.per_file:
                    break
                row = dict(row)
                row["gwas_source"] = tag
                all_rows.append(row)

    if not all_rows:
        print("no rows merged", file=__import__("sys").stderr)
        return 1

    args.outfile.parent.mkdir(parents=True, exist_ok=True)
    with args.outfile.open("w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(all_rows)
    print(f"Wrote {len(all_rows)} rows -> {args.outfile}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
