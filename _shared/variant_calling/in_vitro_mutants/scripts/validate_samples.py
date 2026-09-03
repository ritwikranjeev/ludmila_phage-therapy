#!/usr/bin/env python3
"""
Validate metadata/samples.tsv for the Figure 5a mutation pipeline.
Checks: required columns, file existence, reference_group consistency, no duplicate sample_id.
Exit 0 if valid, 1 otherwise.
"""
import argparse
import csv
import sys
from pathlib import Path


REQUIRED_COLUMNS = [
    "sample_id",
    "sample_type",
    "host_id",
    "reference_group",
    "read1",
    "read2",
]


def main():
    p = argparse.ArgumentParser(description="Validate samples.tsv for Figure 5a pipeline")
    p.add_argument("samples_tsv", help="Path to metadata/samples.tsv")
    p.add_argument("--project-dir", default=".", type=Path, help="Project root for resolving paths")
    p.add_argument("--skip-files", action="store_true", help="Skip checking that read files exist")
    args = p.parse_args()

    project_dir = args.project_dir.resolve()
    samples_path = Path(args.samples_tsv)
    if not samples_path.exists():
        print(f"ERROR: samples file not found: {samples_path}", file=sys.stderr)
        sys.exit(1)

    with open(samples_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        if not rows:
            print("ERROR: samples file is empty or has no data rows", file=sys.stderr)
            sys.exit(1)
        fieldnames = reader.fieldnames or []

    for col in REQUIRED_COLUMNS:
        if col not in fieldnames:
            print(f"ERROR: Missing required column: {col}", file=sys.stderr)
            sys.exit(1)

    errors = []
    sample_ids = [r.get("sample_id", "").strip() for r in rows]
    if len(sample_ids) != len(set(sample_ids)):
        seen = {}
        for i, sid in enumerate(sample_ids):
            if sid in seen:
                errors.append(f"Duplicate sample_id: {sid}")
                break
            seen[sid] = i

    ref_groups = {
        r.get("sample_id", "").strip()
        for r in rows
        if (r.get("sample_type") or "").strip() == "host_reference"
    }
    for r in rows:
        if (r.get("sample_type") or "").strip() != "resistant_mutant":
            continue
        sid = r.get("sample_id", "").strip()
        rg = (r.get("reference_group") or "").strip()
        if rg and rg not in ref_groups:
            errors.append(
                f"resistant_mutant {sid} has reference_group '{rg}' but no host_reference with sample_id '{rg}'"
            )

    if not args.skip_files:
        for r in rows:
            sid = r.get("sample_id", "").strip()
            for key in ("read1", "read2"):
                val = (r.get(key) or "").strip()
                if val in ("", "."):
                    continue
                p = project_dir / val
                if not p.exists():
                    errors.append(f"File not found for {sid} {key}: {p}")

    if errors:
        for e in errors:
            print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    print("Validation passed.", file=sys.stderr)
    sys.exit(0)


if __name__ == "__main__":
    main()
