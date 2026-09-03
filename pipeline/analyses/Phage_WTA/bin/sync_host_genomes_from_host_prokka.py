#!/usr/bin/env python3
"""
Rebuild config/host_genomes.tsv from data/host_prokka/<id>/ directories.

For each host id (subdirectory under data/host_prokka/), finds an assembly FASTA in this
order (first hit wins):
  1. data/hosts/<id>.fa
  2. data/hosts/<id>.fasta
  3. ../../outputs/assemblies_cleaned/<id>/contigs.fasta  (OCPT-Ludmila_Characterization layout)
  4. ../../outputs/assemblies/<id>/contigs.fasta

Paths written to the TSV are relative to the project root when inside the repo tree.

Usage:
  python bin/sync_host_genomes_from_host_prokka.py
  python bin/sync_host_genomes_from_host_prokka.py --dry-run

On HPC3 login nodes with an old system python3, use e.g.:
  module load anaconda/2024.06 && python3 bin/sync_host_genomes_from_host_prokka.py
"""

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Tuple


def project_root() -> Path:
    return Path(__file__).resolve().parent.parent


def assembly_candidates(root: Path, host_id: str) -> List[Path]:
    return [
        root / "data" / "hosts" / f"{host_id}.fa",
        root / "data" / "hosts" / f"{host_id}.fasta",
        root / ".." / ".." / "outputs" / "assemblies_cleaned" / host_id / "contigs.fasta",
        root / ".." / ".." / "outputs" / "assemblies" / host_id / "contigs.fasta",
    ]


def resolve_assembly(root: Path, host_id: str) -> Optional[Path]:
    for p in assembly_candidates(root, host_id):
        try:
            rp = p.resolve()
        except OSError:
            continue
        if rp.is_file():
            return rp
    return None


def rel_or_abs(root: Path, assembly: Path) -> str:
    root_res = root.resolve()
    ar = assembly.resolve()
    try:
        return str(ar.relative_to(root_res)).replace("\\", "/")
    except ValueError:
        return str(ar)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Write config/host_genomes.tsv from data/host_prokka/ directory names.",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="Print id -> path mapping; do not write TSV.",
    )
    args = ap.parse_args()

    root = project_root()
    hp = root / "data" / "host_prokka"
    if not hp.is_dir():
        print(f"ERROR: {hp} is not a directory.", file=sys.stderr)
        return 1

    host_ids = sorted(p.name for p in hp.iterdir() if p.is_dir())
    if not host_ids:
        print(
            f"ERROR: No subdirectories under {hp}. Copy host Prokka outputs first.",
            file=sys.stderr,
        )
        return 1

    rows: List[Tuple[str, str]] = []
    missing: list[str] = []
    for hid in host_ids:
        asm = resolve_assembly(root, hid)
        if asm is None:
            missing.append(hid)
            continue
        rows.append((hid, rel_or_abs(root, asm)))

    if missing:
        print("ERROR: No assembly FASTA found for host id(s):", ", ".join(missing), file=sys.stderr)
        print(
            "Expected one of: data/hosts/<id>.{{fa,fasta}} or "
            "../../outputs/assemblies_cleaned/<id>/contigs.fasta or "
            "../../outputs/assemblies/<id>/contigs.fasta (relative to project root).",
            file=sys.stderr,
        )
        return 1

    out_path = root / "config" / "host_genomes.tsv"
    lines = ["id\tpath"] + [f"{i}\t{p}" for i, p in rows]
    body = "\n".join(lines) + "\n"

    if args.dry_run:
        print(body, end="")
        return 0

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(body, encoding="utf-8")
    print(f"Wrote {out_path} ({len(rows)} host(s)).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
