#!/usr/bin/env python3
"""Post-GWAS: Benjamini–Hochberg FDR on pyseer p-values; sort by q (best first)."""

import argparse
import csv
from pathlib import Path


def bh_fdr(pvals):
    n = len(pvals)
    ranked = sorted(enumerate(pvals), key=lambda x: x[1])
    out = [1.0] * n
    prev = 1.0
    for rank, (idx, p) in reversed(list(enumerate(ranked, start=1))):
        val = min(prev, p * n / rank)
        out[idx] = val
        prev = val
    return out


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Add BH-FDR q-values to pyseer output; sort by q ascending."
    )
    parser.add_argument("--infile", required=True, help="pyseer stdout table")
    parser.add_argument("--outfile", required=True, help="output TSV with qvalue_bh")
    parser.add_argument(
        "--pvalue-col",
        default="lrt-pvalue",
        help="Column to FDR-correct [default: lrt-pvalue]",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=0,
        help="If > 0, keep only the top N rows by q-value after sorting",
    )
    args = parser.parse_args()

    rows = []
    with Path(args.infile).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames
        if not fieldnames:
            print("empty or missing header", file=__import__("sys").stderr)
            return 1
        pcol = args.pvalue_col
        if pcol not in fieldnames:
            # pyseer variants / older builds
            for fallback in ("lrt-pvalue", "filter-pvalue", "pvalue"):
                if fallback in fieldnames:
                    pcol = fallback
                    break
            else:
                print(
                    f"no p-value column (tried {args.pvalue_col}, lrt-pvalue, filter-pvalue, pvalue)",
                    file=__import__("sys").stderr,
                )
                return 1

        for r in reader:
            raw = (r.get(pcol) or "").strip()
            try:
                p = float(raw)
            except ValueError:
                continue
            if p < 0 or p > 1:
                continue
            r["_p_for_fdr"] = p
            rows.append(r)

    if not rows:
        print("no rows with valid p-values", file=__import__("sys").stderr)
        Path(args.outfile).parent.mkdir(parents=True, exist_ok=True)
        Path(args.outfile).write_text("")
        return 0

    qvals = bh_fdr([r["_p_for_fdr"] for r in rows])
    for r, q in zip(rows, qvals):
        r["qvalue_bh"] = q
        del r["_p_for_fdr"]

    rows.sort(key=lambda x: (x["qvalue_bh"], float(x.get(pcol, 1) or 1)))

    if args.top and args.top > 0:
        rows = rows[: args.top]

    extra = ["qvalue_bh"]
    out_fields = [f for f in fieldnames if f is not None] + [
        c for c in extra if c not in (fieldnames or [])
    ]
    Path(args.outfile).parent.mkdir(parents=True, exist_ok=True)
    with Path(args.outfile).open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=out_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows -> {args.outfile}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
