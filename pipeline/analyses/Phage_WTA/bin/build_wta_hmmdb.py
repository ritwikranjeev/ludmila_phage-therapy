#!/usr/bin/env python3
"""
Build per-gene HMMs from config/wta_reference_genes.fasta and a combined DB for hmmsearch/hmmscan.

Usage:
    python bin/build_wta_hmmdb.py

Requires: hmmbuild, hmmpress (HMMER) on PATH.

Modes (see comments in main()):
  - Default: one .hmm per FASTA record under data/wta_hmms/<name>.hmm, then concatenate
    into data/wta_hmms/wta_combined.hmm and run hmmpress (recommended for a single hmmsearch call).
  - Alternative: keep only individual *.hmm files and run hmmsearch once per gene in a driver script
    (more jobs, no hmmpress step).
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path


def project_root() -> Path:
    return Path(__file__).resolve().parent.parent


def read_fasta_records(fasta_path: Path) -> list[tuple[str, str]]:
    """Parse FASTA; skip # comment lines; return (seq_id, sequence) pairs."""
    text = fasta_path.read_text(encoding="utf-8", errors="replace")
    records: list[tuple[str, str]] = []
    current_id: str | None = None
    seq_chunks: list[str] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            if current_id is not None:
                seq = "".join(seq_chunks).replace(" ", "").upper()
                if seq:
                    records.append((current_id, seq))
            current_id = line[1:].strip().split()[0]
            seq_chunks = []
        else:
            seq_chunks.append(line)
    if current_id is not None:
        seq = "".join(seq_chunks).replace(" ", "").upper()
        if seq:
            records.append((current_id, seq))
    return records


def run(cmd: list[str], cwd: Path | None = None) -> None:
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def main() -> int:
    root = project_root()
    ref_fa = root / "config" / "wta_reference_genes.fasta"
    hmm_dir = root / "data" / "wta_hmms"
    combined = hmm_dir / "wta_combined.hmm"

    if not ref_fa.is_file():
        print(f"ERROR: Missing reference FASTA: {ref_fa}", file=sys.stderr)
        return 1

    records = read_fasta_records(ref_fa)
    if not records:
        print("ERROR: No FASTA records found (after skipping # comments).", file=sys.stderr)
        return 1

    hmm_dir.mkdir(parents=True, exist_ok=True)

    # Remove stale combined + press indices so hmmpress can recreate them
    for suffix in ("", ".h3f", ".h3i", ".h3m", ".h3p"):
        p = Path(str(combined) + suffix) if suffix else combined
        if p.is_file():
            p.unlink()

    # Remove prior per-gene HMMs (only known gene stems) to avoid catting old junk
    for seq_id, _ in records:
        old = hmm_dir / f"{seq_id}.hmm"
        if old.is_file():
            old.unlink()

    for seq_id, seq in records:
        if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", seq):
            print(
                f"WARNING: Sequence for {seq_id} contains non-standard AA letters; hmmbuild may fail.",
                file=sys.stderr,
            )
        one_fa = hmm_dir / f".build_{seq_id}.fa"
        one_fa.write_text(f">{seq_id}\n{seq}\n", encoding="utf-8")
        out_hmm = hmm_dir / f"{seq_id}.hmm"
        run(
            [
                "hmmbuild",
                "--amino",
                "-n",
                seq_id,
                str(out_hmm),
                str(one_fa),
            ]
        )
        one_fa.unlink(missing_ok=True)

    # --- Combined DB (optional mode B: comment out this block and use individual *.hmm only) ---
    parts = [hmm_dir / f"{sid}.hmm" for sid, _ in records]
    for p in parts:
        if not p.is_file():
            print(f"ERROR: Expected HMM missing: {p}", file=sys.stderr)
            return 1
    with combined.open("wb") as out:
        for p in parts:
            out.write(p.read_bytes())
    run(["hmmpress", "-f", str(combined.name)], cwd=hmm_dir)
    # hmmpress -f overwrites; leaves wta_combined.hmm + .h3*

    print(f"Done. Individual HMMs under {hmm_dir}/<gene>.hmm")
    print(f"Combined + pressed DB: {combined} (+ .h3*)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
