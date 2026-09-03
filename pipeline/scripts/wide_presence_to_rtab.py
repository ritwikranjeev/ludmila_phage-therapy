#!/usr/bin/env python3
"""Transpose wide sample×feature TSV to pyseer/roary-style Rtab (Gene × samples, 0/1)."""

import csv
import sys


def main():
    if len(sys.argv) != 2:
        print("usage: wide_presence_to_rtab.py <in.tsv> > out.Rtab", file=sys.stderr)
        return 1
    with open(sys.argv[1], newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        print("empty input", file=sys.stderr)
        return 1
    if "sample" not in rows[0]:
        print("need a 'sample' column in header", file=sys.stderr)
        return 1
    samples = [r["sample"].strip() for r in rows]
    features = [k for k in rows[0].keys() if k != "sample"]
    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(["Gene"] + samples)
    for feat in features:
        out = [feat]
        for r in rows:
            v = (r.get(feat, "0") or "0").strip()
            if v not in ("0", "1"):
                v = "1" if v and v not in ("0", "false", "False") else "0"
            out.append(v)
        writer.writerow(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
