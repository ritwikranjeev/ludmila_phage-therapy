#!/usr/bin/env python3
"""
One-shot local troubleshooting gate for Phage_WTA (no Slurm).

Run from project root:
    python3 bin/validate_troubleshoot_complete.py

Exit 0: syntax OK, Slurm→bin wiring OK, smoke-equivalent Python steps passed.
Exit 1: any check failed (see stderr).

Optional:
    python3 bin/validate_troubleshoot_complete.py --print-hpc3-commands
    # prints a minimal SSH checklist (you run on your Mac terminal).

    Override defaults with env vars:
      PHAGE_WTA_HPC3_USER, PHAGE_WTA_HPC3_HOST, PHAGE_WTA_HPC3_PROJECT_DIR
"""

from __future__ import annotations

import argparse
import ast
import os
import re
import subprocess
import sys
from pathlib import Path


def project_root() -> Path:
    return Path(__file__).resolve().parent.parent


def check_python_syntax(bin_dir: Path) -> list[str]:
    errors: list[str] = []
    for py in sorted(bin_dir.glob("*.py")):
        try:
            ast.parse(py.read_text(encoding="utf-8", errors="replace"), filename=str(py))
        except SyntaxError as e:
            errors.append(f"SyntaxError {py}: {e}")
    return errors


def check_bash_syntax(root: Path) -> list[str]:
    errors: list[str] = []
    for folder in ("slurm", "bin"):
        d = root / folder
        if not d.is_dir():
            errors.append(f"Missing directory {d}")
            continue
        for sh in sorted(d.glob("*.sh")):
            r = subprocess.run(
                ["bash", "-n", str(sh)],
                capture_output=True,
                text=True,
            )
            if r.returncode != 0:
                errors.append(f"bash -n {sh}: {r.stderr or r.stdout or 'failed'}")
    return errors


def slurm_referenced_python_scripts(root: Path) -> list[str]:
    """Return basenames of bin/*.py invoked from slurm/*.sh (real invocations only)."""
    slurm = root / "slurm"
    found: set[str] = set()
    # python3 "${ROOT}/bin/foo.py" or python3 "${PHAGE_WTA_PROJECT_ROOT}/bin/foo.py"
    pat = re.compile(
        r'python3\s+["\']?\$\{PHAGE_WTA_PROJECT_ROOT\}/bin/([A-Za-z0-9_]+\.py)["\']?'
    )
    for sh in sorted(slurm.glob("*.sh")):
        text = sh.read_text(encoding="utf-8", errors="replace")
        for m in pat.finditer(text):
            found.add(m.group(1))
    return sorted(found)


def check_slurm_bin_wiring(root: Path) -> list[str]:
    errors: list[str] = []
    for name in slurm_referenced_python_scripts(root):
        p = root / "bin" / name
        if not p.is_file():
            errors.append(f"Slurm references missing {p}")
    return errors


def run_smoke_steps(root: Path) -> list[str]:
    errors: list[str] = []
    env = {**os.environ, "PHAGE_WTA_PROJECT_ROOT": str(root), "PHAGE_WTA_SKIP_CONDA": "1"}

    def run(args: list[str]) -> None:
        r = subprocess.run(args, cwd=str(root), env=env, capture_output=True, text=True)
        if r.returncode != 0:
            msg = (r.stderr or r.stdout or "").strip() or "(no output)"
            errors.append(f"Command {' '.join(args)} exit {r.returncode}:\n{msg}")

    run(["bash", str(root / "bin" / "sanity_check_paths.sh")])
    run([sys.executable, str(root / "bin" / "prepare_genome_lists.py")])
    r = subprocess.run(
        [sys.executable, str(root / "bin" / "preflight_check.py"), "--allow-wta-placeholder"],
        cwd=str(root),
        env=env,
        capture_output=True,
        text=True,
    )
    if r.returncode != 0:
        msg = (r.stderr or r.stdout or "").strip() or "(no output)"
        errors.append(f"preflight_check.py --allow-wta-placeholder exit {r.returncode}:\n{msg}")
    run(
        [
            sys.executable,
            str(root / "bin" / "build_wta_feature_matrix.py"),
            "--allow-empty-raw",
        ]
    )
    return errors


def print_hpc3_commands() -> None:
    user = os.environ.get("PHAGE_WTA_HPC3_USER", "ritwik1")
    host = os.environ.get("PHAGE_WTA_HPC3_HOST", "hpc3.rcic.uci.edu")
    base = os.environ.get(
        "PHAGE_WTA_HPC3_PROJECT_DIR",
        "/pub/ritwik1/manuscripts/OCPT-Ludmila_PhageTherapy/WTA-RBP_Analysis",
    )
    print(
        f"""# Run on your Mac (SSH + project on HPC3).
# Defaults: user={user!r}, host={host!r}, cd={base!r}
# Override: export PHAGE_WTA_HPC3_USER PHAGE_WTA_HPC3_HOST PHAGE_WTA_HPC3_PROJECT_DIR

ssh -o BatchMode=yes -o ConnectTimeout=25 {user}@{host} bash -lc '
set -euo pipefail
cd "{base}"
bash bin/sanity_check_paths.sh
python3 bin/preflight_check.py --allow-wta-placeholder
export DRY_RUN=1
jid_p=$(sbatch --parsable --array=1-1 slurm/annotate_phage_prokka.sh)
jid_h=$(sbatch --parsable --array=1-1 slurm/annotate_hosts_prokka.sh)
echo "Submitted phage job $jid_p host job $jid_h — wait then: ls -la logs/phage_annotate_*_1.out logs/host_annotate_*_1.out"
'
"""
    )


def main() -> int:
    ap = argparse.ArgumentParser(description="Phage_WTA local troubleshoot gate.")
    ap.add_argument(
        "--print-hpc3-commands",
        action="store_true",
        help=(
            "Print SSH/sbatch checklist for HPC3 (does not run checks). "
            "Env: PHAGE_WTA_HPC3_USER, PHAGE_WTA_HPC3_HOST, PHAGE_WTA_HPC3_PROJECT_DIR."
        ),
    )
    args = ap.parse_args()
    if args.print_hpc3_commands:
        print_hpc3_commands()
        return 0

    root = project_root()
    all_err: list[str] = []

    all_err.extend(check_python_syntax(root / "bin"))
    all_err.extend(check_bash_syntax(root))
    all_err.extend(check_slurm_bin_wiring(root))
    all_err.extend(run_smoke_steps(root))

    if all_err:
        print("validate_troubleshoot_complete.py FAILED:", file=sys.stderr)
        for e in all_err:
            print(e, file=sys.stderr)
        return 1

    print("validate_troubleshoot_complete.py PASSED (local gate).")
    print("HPC3: run with --print-hpc3-commands and execute that block when ready.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
