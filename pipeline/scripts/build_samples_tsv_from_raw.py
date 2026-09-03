#!/usr/bin/env python3
"""Emit config/samples.tsv by scanning data/raw/ on disk.

Rules:
  - Exactly 2 FASTQs -> platform illumina, r1/r2 full paths (R1/R2 or _1/_2).
  - Exactly 1 FASTQ  -> platform long, longread path, r1/r2 = NA.

Run on the machine that holds data/raw (e.g. HPC), from project root:
  python3 scripts/build_samples_tsv_from_raw.py > config/samples.tsv

Optional:
  RAW_ROOT=/other/raw python3 scripts/build_samples_tsv_from_raw.py
"""

import os
import sys
from pathlib import Path
from typing import List, Tuple

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_ROOT = Path(os.environ.get("RAW_ROOT", PROJECT_ROOT / "data" / "raw"))

def list_fastqs(d: Path) -> List[Path]:
    if not d.is_dir():
        return []
    out = []
    for p in d.iterdir():
        if not p.is_file():
            continue
        n = p.name.lower()
        if n.endswith((".fastq.gz", ".fastq", ".fq.gz", ".fq")):
            out.append(p)
    return sorted(out, key=lambda x: x.name)


def is_r1(name: str) -> bool:
    n = name.lower()
    return "_r1" in n or "_1.fastq" in n or "_1.fq" in n


def is_r2(name: str) -> bool:
    n = name.lower()
    return "_r2" in n or "_2.fastq" in n or "_2.fq" in n


def order_pe_pair(a: Path, b: Path) -> Tuple[Path, Path]:
    an, bn = a.name, b.name
    if is_r1(an) and is_r2(bn):
        return a, b
    if is_r1(bn) and is_r2(an):
        return b, a
    return tuple(sorted((a, b), key=lambda x: x.name))


def norm_path(p: Path) -> str:
    """Prefer /pub/... over /dfs6b/pub/... for docs and shared paths on HPC3."""
    s = str(p.resolve())
    if "/dfs6b/pub/" in s:
        s = s.replace("/dfs6b/pub/", "/pub/", 1)
    return s


def main() -> None:
    if not RAW_ROOT.is_dir():
        print(f"ERROR: RAW_ROOT not found: {RAW_ROOT}", file=sys.stderr)
        sys.exit(1)

    rows = []  # type: List[Tuple[str, str, str, str, str]]
    for sample_dir in sorted(RAW_ROOT.iterdir(), key=lambda x: x.name):
        if not sample_dir.is_dir():
            continue
        sample = sample_dir.name
        fqs = list_fastqs(sample_dir)
        if len(fqs) == 1:
            rows.append((sample, "long", "NA", "NA", norm_path(fqs[0])))
        elif len(fqs) == 2:
            r1, r2 = order_pe_pair(fqs[0], fqs[1])
            rows.append(
                (sample, "illumina", norm_path(r1), norm_path(r2), "NA")
            )
        else:
            print(f"WARN skip {sample}: {len(fqs)} fastq files", file=sys.stderr)

    print("sample\tplatform\tr1\tr2\tlongread")
    for r in rows:
        print("\t".join(r))


if __name__ == "__main__":
    main()
