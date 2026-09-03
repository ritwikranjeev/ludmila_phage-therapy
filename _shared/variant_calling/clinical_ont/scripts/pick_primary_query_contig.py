#!/usr/bin/env python3
"""Pick the post (query) contig with the most reference-aligned bases to a given ref contig (show-coords)."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path
from collections import defaultdict


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("delta", help="nucmer .delta file")
    p.add_argument("ref_contig", help="Reference sequence id (e.g. contig_2)")
    args = p.parse_args(argv)

    sc = "show-coords"
    prefix = os.environ.get("CONDA_PREFIX", "")
    if prefix:
        cand = Path(prefix) / "bin" / "show-coords"
        if cand.is_file():
            sc = str(cand)
    try:
        proc = subprocess.run(
            [sc, "-H", "-c", "-l", args.delta],
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError:
        print("show-coords not on PATH (need MUMmer4 conda env)", file=sys.stderr)
        return 1
    except subprocess.CalledProcessError as e:
        print(e.stderr, file=sys.stderr)
        return e.returncode

    # Per alignment row, use max of the third | block (alignment length in coords);
    # then pick query contig with the largest single block (main syntenic segment).
    best_block: dict[str, int] = defaultdict(int)
    for line in proc.stdout.splitlines():
        line = line.strip()
        if not line or line[0] not in "0123456789":
            continue
        parts = [x.strip() for x in line.split("|")]
        if len(parts) < 7:
            continue
        tail = parts[-1].split()
        if len(tail) < 2:
            continue
        ref, qry = tail[-2], tail[-1]
        if ref != args.ref_contig:
            continue
        try:
            tnums = [int(x) for x in parts[2].split()]
            if len(tnums) < 2:
                continue
            block = max(tnums[0], tnums[1])
        except (IndexError, ValueError):
            continue
        if block > best_block[qry]:
            best_block[qry] = block

    if not best_block:
        print(f"No alignments to {args.ref_contig!r} in delta", file=sys.stderr)
        return 1
    best = max(best_block.items(), key=lambda kv: kv[1])
    print(best[0])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
