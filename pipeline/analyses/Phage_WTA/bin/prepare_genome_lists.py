#!/usr/bin/env python3
"""Scan data/phage and data/hosts for FASTA files; write config/phage_genomes.tsv and config/host_genomes.tsv.

Paths under the project root are stored relative to that root (POSIX slashes) so lists commit cleanly;
Slurm drivers resolve them against PHAGE_WTA_PROJECT_ROOT. Paths outside the tree stay absolute.
"""

from __future__ import annotations

import sys
from pathlib import Path

FASTA_SUFFIXES = {".fasta", ".fa"}


def project_root() -> Path:
    return Path(__file__).resolve().parent.parent


def list_fasta_files(directory: Path, label: str) -> list[Path]:
    if not directory.exists():
        print(
            f"WARNING: {label} directory does not exist: {directory}",
            file=sys.stderr,
        )
        return []
    if not directory.is_dir():
        print(
            f"WARNING: {label} path is not a directory: {directory}",
            file=sys.stderr,
        )
        return []

    files: list[Path] = []
    for path in sorted(directory.iterdir()):
        if not path.is_file():
            continue
        if path.suffix.lower() not in FASTA_SUFFIXES:
            continue
        files.append(path)
    return files


def write_tsv(out_path: Path, rows: list[tuple[str, str]]) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as fh:
        fh.write("id\tpath\n")
        for genome_id, abs_path in rows:
            fh.write(f"{genome_id}\t{abs_path}\n")


def build_rows(root: Path, files: list[Path], label: str) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    seen: set[str] = set()
    root_res = root.resolve()
    for path in files:
        genome_id = path.stem
        if genome_id in seen:
            print(
                f"WARNING: duplicate id '{genome_id}' in {label}, skipping {path}",
                file=sys.stderr,
            )
            continue
        seen.add(genome_id)
        abs_p = path.resolve()
        try:
            rel = abs_p.relative_to(root_res)
            path_str = str(rel).replace("\\", "/")
        except ValueError:
            path_str = str(abs_p)
        rows.append((genome_id, path_str))
    return rows


def main() -> int:
    root = project_root()
    phage_dir = root / "data" / "phage"
    hosts_dir = root / "data" / "hosts"
    config_dir = root / "config"

    phage_files = list_fasta_files(phage_dir, "Phage genomes (data/phage)")
    host_files = list_fasta_files(hosts_dir, "Host genomes (data/hosts)")

    phage_rows = build_rows(root, phage_files, "data/phage")
    host_rows = build_rows(root, host_files, "data/hosts")

    write_tsv(config_dir / "phage_genomes.tsv", phage_rows)
    write_tsv(config_dir / "host_genomes.tsv", host_rows)

    print(f"Wrote {len(phage_rows)} phage row(s) to {config_dir / 'phage_genomes.tsv'}")
    print(f"Wrote {len(host_rows)} host row(s) to {config_dir / 'host_genomes.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
