#!/usr/bin/env python3
"""Summarize per-sample prophage output tables into one TSV."""

import csv
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
IN_ROOT = PROJECT_ROOT / "outputs" / "prophages"
OUT_TSV = PROJECT_ROOT / "outputs" / "prophages" / "prophage_summary.tsv"
FIELDS = [
    "sample", "prophage_id", "contig", "start", "end", "size_bp",
    "score", "closest_known_phage",
]


def _append_from_prophage_coordinates(rows, sample, path):
    """PhiSpy writes prophage_coordinates.tsv with no header: id, contig, start, stop, att…"""
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            pid, contig, start_s, end_s = parts[0], parts[1], parts[2], parts[3]
            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                continue
            lo, hi = (start, end) if start <= end else (end, start)
            rows.append(
                {
                    "sample": sample,
                    "prophage_id": pid,
                    "contig": contig,
                    "start": lo,
                    "end": hi,
                    "size_bp": max(0, hi - lo + 1),
                    "score": "",
                    "closest_known_phage": "",
                }
            )


def _append_from_prophage_tsv(rows, sample, path):
    with path.open() as handle:
        rdr = csv.DictReader(handle, delimiter="\t")
        for r in rdr:
            start = int(r.get("Start", r.get("start", "0")))
            end = int(r.get("Stop", r.get("end", "0")))
            lo, hi = (start, end) if start <= end else (end, start)
            pid = r.get("Prophage number", r.get("prophage_id", ""))
            contig = r.get("Contig", r.get("contig", ""))
            rows.append(
                {
                    "sample": sample,
                    "prophage_id": pid,
                    "contig": contig,
                    "start": lo,
                    "end": hi,
                    "size_bp": max(0, hi - lo + 1),
                    "score": r.get("score", ""),
                    "closest_known_phage": r.get("closest_known_phage", ""),
                }
            )


def main():
    rows = []
    for sample_dir in sorted(IN_ROOT.glob("*")):
        if not sample_dir.is_dir():
            continue
        sample = sample_dir.name
        coords = sample_dir / "prophage_coordinates.tsv"
        simple = sample_dir / "prophage.tsv"
        legacy = sample_dir / "prophage_calls.tsv"
        if coords.exists():
            _append_from_prophage_coordinates(rows, sample, coords)
        elif simple.exists():
            _append_from_prophage_tsv(rows, sample, simple)
        elif legacy.exists():
            with legacy.open() as handle:
                for r in csv.DictReader(handle, delimiter="\t"):
                    start = int(r.get("start", "0"))
                    end = int(r.get("end", "0"))
                    lo, hi = (start, end) if start <= end else (end, start)
                    rows.append(
                        {
                            "sample": sample,
                            "prophage_id": r.get("prophage_id", ""),
                            "contig": r.get("contig", ""),
                            "start": lo,
                            "end": hi,
                            "size_bp": max(0, hi - lo + 1),
                            "score": r.get("score", ""),
                            "closest_known_phage": r.get("closest_known_phage", ""),
                        }
                    )
    with OUT_TSV.open("w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {OUT_TSV} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
