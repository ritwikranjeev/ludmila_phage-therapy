#!/usr/bin/env python3
"""
Verify Phage_WTA project layout and toolchain before HPC3 runs.

Does not source conda; uses the current interpreter and PATH.

Exit 0: critical checks pass (some optional tools may be missing — see summary).
Exit 1: critical failure (missing config, genome TSVs, python issues, or placeholder WTA refs).
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

if sys.version_info < (3, 9):
    print(
        "ERROR: preflight_check.py needs Python >= 3.9.\n"
        "On UCI HPC3, login-node `python3` is often an older user-managed interpreter.\n"
        "Use: bash bin/hpc3_login_preflight.sh [...]   or: module load anaconda/2024.06 && python3 ...",
        file=sys.stderr,
    )
    raise SystemExit(1)

CRITICAL_DIRS = [
    "logs",
    "data/phage",
    "data/hosts",
    "data/phage_prokka",
    "data/host_prokka",
    "data/phage_rbp",
    "data/wta_features",
    "data/wta_features/by_host",
    "data/wta_hmms",
    "data/alphafold/Ludmila_RBPs",
    "data/phenotypes",
    "analysis",
]

TOOLS_OPTIONAL = [
    ("prokka", "Prokka"),
    ("hmmbuild", "HMMER hmmbuild"),
    ("hmmsearch", "HMMER hmmsearch"),
    ("hmmpress", "HMMER hmmpress"),
    ("blastp", "BLAST+ blastp"),
    ("colabfold_batch", "ColabFold batch"),
]


def project_root() -> Path:
    return Path(__file__).resolve().parent.parent


def read_genome_tsv_ids(tsv_path: Path) -> list[str]:
    text = tsv_path.read_text(encoding="utf-8", errors="replace").strip()
    if not text:
        return []
    ids: list[str] = []
    for line in text.splitlines()[1:]:
        line = line.strip()
        if not line:
            continue
        hid = line.split("\t", 1)[0].strip()
        if hid:
            ids.append(hid)
    return ids


def host_prokka_subdir_ids(root: Path) -> set[str]:
    hp = root / "data" / "host_prokka"
    if not hp.is_dir():
        return set()
    return {p.name for p in hp.iterdir() if p.is_dir()}


def phage_prokka_subdir_ids(root: Path) -> set[str]:
    pp = root / "data" / "phage_prokka"
    if not pp.is_dir():
        return set()
    return {p.name for p in pp.iterdir() if p.is_dir()}


def warn_annotation_layout_mismatch(root: Path, warn: list[str]) -> None:
    """If annotation folders exist, flag TSV ids that do not match (common HPC3 mistake)."""
    hg = root / "config" / "host_genomes.tsv"
    pg = root / "config" / "phage_genomes.tsv"
    hp_ids = host_prokka_subdir_ids(root)
    tsv_host = set(read_genome_tsv_ids(hg)) if hg.is_file() else set()
    if hp_ids:
        only_dirs = sorted(hp_ids - tsv_host)
        only_tsv = sorted(tsv_host - hp_ids)
        if only_dirs or only_tsv:
            if only_dirs:
                warn.append(
                    f"host_prokka/ has {len(only_dirs)} id(s) missing from host_genomes.tsv "
                    f"(e.g. {', '.join(only_dirs[:5])}{'…' if len(only_dirs) > 5 else ''}); "
                    f"run: python3 bin/sync_host_genomes_from_host_prokka.py"
                )
            if only_tsv:
                warn.append(
                    f"host_genomes.tsv lists id(s) with no data/host_prokka/<id>/: "
                    f"{', '.join(only_tsv[:8])}{'…' if len(only_tsv) > 8 else ''}"
                )

    pp_ids = phage_prokka_subdir_ids(root)
    tsv_phage = set(read_genome_tsv_ids(pg)) if pg.is_file() else set()
    if pp_ids:
        only_dirs = sorted(pp_ids - tsv_phage)
        only_tsv = sorted(tsv_phage - pp_ids)
        if only_dirs or only_tsv:
            if only_dirs:
                warn.append(
                    f"phage_prokka/ id(s) not in phage_genomes.tsv: {', '.join(only_dirs)}"
                )
            if only_tsv:
                warn.append(
                    f"phage_genomes.tsv id(s) with no data/phage_prokka/<id>/: "
                    f"{', '.join(only_tsv)}"
                )


def main() -> int:
    ap = argparse.ArgumentParser(description="Phage_WTA preflight checks.")
    ap.add_argument(
        "--allow-wta-placeholder",
        action="store_true",
        help="Do not fail if MTESTSEQ remains in wta_reference_genes.fasta (smoke tests only).",
    )
    args = ap.parse_args()

    root = project_root()
    fail: list[str] = []
    warn: list[str] = []

    env_setup = root / "config" / "env_setup.sh"
    if not env_setup.is_file():
        fail.append(f"Missing {env_setup}")
    else:
        print(f"OK  config/env_setup.sh")

    for name in ("config/phage_genomes.tsv", "config/host_genomes.tsv"):
        p = root / name
        if not p.is_file():
            fail.append(f"Missing {p}")
        else:
            text = p.read_text(encoding="utf-8", errors="replace").strip()
            if not text or "\t" not in text.splitlines()[0]:
                fail.append(f"Invalid or empty {p}")
            else:
                data_lines = [ln for ln in text.splitlines()[1:] if ln.strip()]
                if not data_lines:
                    fail.append(f"No data rows in {p}")
                else:
                    print(f"OK  {name} ({len(data_lines)} data row(s))")

    wta_fa = root / "config" / "wta_reference_genes.fasta"
    if not wta_fa.is_file():
        fail.append(f"Missing {wta_fa}")
    else:
        body = wta_fa.read_text(encoding="utf-8", errors="replace")
        if "MTESTSEQ" in body.upper():
            msg = (
                "config/wta_reference_genes.fasta still contains MTESTSEQ placeholder; "
                "replace with real reference protein sequences before HMM builds."
            )
            if args.allow_wta_placeholder:
                warn.append(msg + " (--allow-wta-placeholder)")
                print("WARN WTA reference still has MTESTSEQ (allowed by flag)")
            else:
                fail.append(msg)
        else:
            print("OK  WTA reference FASTA (no MTESTSEQ placeholder detected)")

    py = shutil.which("python3") or shutil.which("python")
    if not py:
        fail.append("python3/python not found on PATH")
    else:
        print(f"OK  interpreter: {py}")

    for relp in CRITICAL_DIRS:
        p = root / relp
        try:
            p.mkdir(parents=True, exist_ok=True)
        except OSError as e:
            fail.append(f"Cannot create {p}: {e}")
        if p.is_dir():
            print(f"OK  dir {relp}/")
        else:
            fail.append(f"Not a directory after mkdir: {p}")

    for exe, label in TOOLS_OPTIONAL:
        path = shutil.which(exe)
        if path:
            print(f"OK  {label}: {path}")
        else:
            warn.append(f"MISSING (PATH) {label} (`{exe}`)")

    warn_annotation_layout_mismatch(root, warn)

    print("\n--- Preflight summary ---")
    if warn:
        print("Warnings (install or module-load before real runs):")
        for w in warn:
            print(f"  {w}")
    if fail:
        print("CRITICAL failures:")
        for f in fail:
            print(f"  {f}")
        print("\nPreflight FAILED.")
        return 1

    print("Preflight PASSED (critical checks).")
    if warn:
        print("Resolve warnings above for full pipeline execution.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
