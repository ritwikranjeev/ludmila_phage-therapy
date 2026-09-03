#!/usr/bin/env python3
"""
Create or validate data/phenotypes/ludmila_infectivity.tsv for WTA–infectivity association.

Columns: host_id, infectivity (0/1 or numeric efficiency; downstream script binarizes for logistic regression).

Lines starting with # are comments (ignored by pandas in associate_wta_with_infectivity.py).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REQUIRED_COLS = ("host_id", "infectivity")

EXAMPLE_LINES = """# DUMMY DATA — Replace every row below with real Ludmila infectivity measurements.
# host_id values must match host_id in data/wta_features/wta_feature_matrix.tsv.
host_id	infectivity
strain_A	1
strain_B	0
strain_C	1
strain_D	0.75
"""


def project_root() -> Path:
    return Path(__file__).resolve().parent.parent


def write_from_matrix(matrix_path: Path, out_path: Path) -> int:
    """Write host_id + empty infectivity column from wta_feature_matrix.tsv header column host_id."""
    text = matrix_path.read_text(encoding="utf-8", errors="replace")
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        print(f"ERROR: Empty matrix: {matrix_path}", file=sys.stderr)
        return 1
    header = lines[0].split("\t")
    if "host_id" not in header:
        print(f"ERROR: No host_id column in {matrix_path}", file=sys.stderr)
        return 1
    col = header.index("host_id")
    host_ids: list[str] = []
    for line in lines[1:]:
        parts = line.split("\t")
        if col < len(parts) and parts[col].strip():
            host_ids.append(parts[col].strip())
    if not host_ids:
        print(f"ERROR: No host rows in {matrix_path}", file=sys.stderr)
        return 1
    body = "\n".join(
        [
            "# Fill infectivity with 0 or 1 per host (numeric OK; association binarizes >0).",
            "# Rows generated from wta_feature_matrix.tsv — keep host_id spelling identical.",
            "# Then replace data/phenotypes/ludmila_infectivity.tsv or point your workflow at this file.",
            "host_id\tinfectivity",
            *[f"{hid}\t" for hid in host_ids],
        ]
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(body + "\n", encoding="utf-8")
    print(f"Wrote {out_path} ({len(host_ids)} host row(s), infectivity column empty — fill before association).")
    return 0


def validate(path: Path) -> int:
    text = path.read_text(encoding="utf-8")
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        print(f"ERROR: No data rows in {path}", file=sys.stderr)
        return 1
    header = lines[0].split("\t")
    for col in REQUIRED_COLS:
        if col not in header:
            print(f"ERROR: Missing column {col!r} in {path}", file=sys.stderr)
            return 1
    print(f"OK: {path} validates (header includes {REQUIRED_COLS}).")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="Create or check ludmila_infectivity.tsv.")
    ap.add_argument(
        "--path",
        type=Path,
        default=None,
        help="TSV path (default: data/phenotypes/ludmila_infectivity.tsv under project root).",
    )
    ap.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing file with the example template (use with care).",
    )
    ap.add_argument(
        "--from-matrix",
        action="store_true",
        help="Write a TSV with host_id from data/wta_features/wta_feature_matrix.tsv and empty infectivity.",
    )
    ap.add_argument(
        "--matrix-path",
        type=Path,
        default=None,
        help="Matrix TSV for --from-matrix (default: data/wta_features/wta_feature_matrix.tsv).",
    )
    ap.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output path for --from-matrix (default: data/phenotypes/ludmila_infectivity.TO_FILL.tsv).",
    )
    args = ap.parse_args()

    root = project_root()

    if args.from_matrix:
        matrix = (
            args.matrix_path
            if args.matrix_path is not None
            else root / "data" / "wta_features" / "wta_feature_matrix.tsv"
        ).resolve()
        if not matrix.is_file():
            print(f"ERROR: Matrix not found: {matrix}", file=sys.stderr)
            return 1
        out = (
            args.output
            if args.output is not None
            else root / "data" / "phenotypes" / "ludmila_infectivity.TO_FILL.tsv"
        ).resolve()
        if out.is_file() and not args.force:
            print(f"ERROR: {out} exists; use --force to overwrite.", file=sys.stderr)
            return 1
        return write_from_matrix(matrix, out)

    path = (args.path if args.path is not None else root / "data" / "phenotypes" / "ludmila_infectivity.tsv").resolve()

    if path.is_file() and not args.force:
        return validate(path)

    if path.is_file() and args.force:
        path.unlink()

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(EXAMPLE_LINES, encoding="utf-8")
    print(f"Wrote example template to {path}")
    print("Replace dummy rows with real Ludmila infectivity data before association analysis.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
