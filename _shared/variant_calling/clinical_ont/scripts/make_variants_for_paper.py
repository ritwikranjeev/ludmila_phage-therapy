#!/usr/bin/env python3
"""Compose manuscript-ready variant table from summary + WTA/envelope subset."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--summary", type=Path, required=True)
    p.add_argument("--focus", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()

    summary = read_tsv(args.summary)
    focus = read_tsv(args.focus)
    focus_keys = {(r["chrom"], r["pos"], r["ref"], r["alt"]) for r in focus}

    keep = []
    for row in summary:
        key = (row["chrom"], row["pos"], row["ref"], row["alt"])
        tier = row.get("tier", "")
        if tier not in {"Tier 1", "Tier 2"}:
            continue
        row["focus_set"] = "yes" if key in focus_keys else "no"
        keep.append(row)

    keep.sort(
        key=lambda r: (
            r.get("tier", ""),
            r.get("category", ""),
            r.get("chrom", ""),
            int(r.get("pos", "0") or 0),
        )
    )

    fields = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "gene",
        "locus_tag",
        "product",
        "effect",
        "impact",
        "category",
        "tier",
        "focus_set",
        "caller_support",
        "pass_both_callers",
        "medaka_qual",
        "medaka_dp",
        "medaka_af",
        "clair3_qual",
        "clair3_dp",
        "clair3_af",
        "region_flag",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as oh:
        w = csv.DictWriter(oh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for row in keep:
            w.writerow({k: row.get(k, "") for k in fields})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
