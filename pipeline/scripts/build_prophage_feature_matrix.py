#!/usr/bin/env python3
"""Build sample×binary-feature TSV for pyseer from per-sample PhiSpy coordinates."""

import csv
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SAMPLES_TSV = PROJECT_ROOT / "config" / "samples.tsv"
PROPHAGE_ROOT = PROJECT_ROOT / "outputs" / "prophages"
OUT_PATH = PROJECT_ROOT / "outputs" / "prophages" / "prophage_feature_matrix.tsv"

# At least one region this long (bp) counts as a "large" prophage for a binary feature.
LARGE_MIN_BP = 25000


def parse_coordinates(path: Path):
    """Return (n_regions, max_size_bp) from prophage_coordinates.tsv (no header)."""
    n = 0
    max_sz = 0
    if not path.is_file():
        return 0, 0
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            try:
                a = int(parts[2])
                b = int(parts[3])
            except ValueError:
                continue
            lo, hi = (a, b) if a <= b else (b, a)
            max_sz = max(max_sz, hi - lo + 1)
            n += 1
    return n, max_sz


def main():
    samples = []
    with SAMPLES_TSV.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            s = (row.get("sample") or "").strip()
            if s:
                samples.append(s)

    fields = ["sample", "pp_any", "pp_multi", "pp_large_any"]
    rows_out = []
    for sample in samples:
        coords = PROPHAGE_ROOT / sample / "prophage_coordinates.tsv"
        n, max_sz = parse_coordinates(coords)
        rows_out.append(
            {
                "sample": sample,
                "pp_any": 1 if n >= 1 else 0,
                "pp_multi": 1 if n >= 2 else 0,
                "pp_large_any": 1 if max_sz >= LARGE_MIN_BP else 0,
            }
        )

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUT_PATH.open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows_out)
    print(f"Wrote {OUT_PATH} ({len(rows_out)} samples, features: {fields[1:]})")


if __name__ == "__main__":
    main()
