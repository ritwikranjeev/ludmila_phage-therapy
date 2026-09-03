#!/usr/bin/env python3
"""Report missing output.gd for resistant_mutants; optional --rm-stale removes incomplete breseq trees."""
import argparse
import csv
import shutil
import sys
from pathlib import Path


def mutant_ids(samples_tsv):
    with open(samples_tsv, newline="") as f:
        return [
            r["sample_id"].strip()
            for r in csv.DictReader(f, delimiter="\t")
            if (r.get("sample_type") or "").strip() == "resistant_mutant"
        ]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("samples_tsv", help="metadata/samples.tsv")
    p.add_argument(
        "--results-dir",
        type=Path,
        required=True,
        help="Path to results/ (contains breseq/)",
    )
    p.add_argument(
        "--rm-stale",
        action="store_true",
        help="Remove results/breseq/<id> when output/output.gd is missing but directory exists",
    )
    args = p.parse_args()

    ids = mutant_ids(args.samples_tsv)
    br = args.results_dir / "breseq"

    def scan():
        missing, stale, ok = [], [], []
        for sid in ids:
            gd = br / sid / "output" / "output.gd"
            d = br / sid
            if gd.is_file():
                ok.append(sid)
            elif d.is_dir():
                stale.append(sid)
            else:
                missing.append(sid)
        return missing, stale, ok

    missing, stale, ok = scan()
    if args.rm_stale and stale:
        for sid in stale:
            d = br / sid
            shutil.rmtree(d, ignore_errors=True)
            print(f"removed stale {d}", file=sys.stderr)
        missing, stale, ok = scan()

    print(f"OK (have output.gd): {len(ok)} / {len(ids)}")
    if missing:
        print(f"Not started (no dir): {len(missing)} -> {', '.join(missing)}")
    if stale:
        print(f"Stale (dir but no .gd): {len(stale)} -> {', '.join(stale)}")

    return 0 if len(ok) == len(ids) else 1


if __name__ == "__main__":
    sys.exit(main())
