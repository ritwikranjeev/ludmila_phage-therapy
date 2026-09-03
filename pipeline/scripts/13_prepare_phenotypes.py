#!/usr/bin/env python3
"""Prepare clean phenotype table from host_range.tsv.

Derived fields:
  strong_infect, weak_infect, infection_class — primary phenotypes for association.
  group — optional lineage/metadata label (BOR, DTD, VISA, USA300, other) for
          descriptive plots or covariates only; not a primary analysis stratum.
"""

import csv
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SAMPLES_TSV = PROJECT_ROOT / "config" / "samples.tsv"
IN_PHENO = PROJECT_ROOT / "data" / "phenotypes" / "host_range.tsv"
OUT_PHENO = PROJECT_ROOT / "config" / "phenotypes.tsv"


def infer_group(sample):
    """Optional label from sample prefix (covariates / faceting only)."""
    prefix = sample.split("-")[0]
    return prefix if prefix in {"BOR", "DTD", "VISA", "USA300"} else "other"


def main():
    samples = set()
    with SAMPLES_TSV.open() as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            samples.add(r["sample"])

    out_rows = []
    with IN_PHENO.open() as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            sample = r["sample"]
            if sample not in samples:
                raise ValueError(f"Phenotype sample not in samples.tsv: {sample}")
            d = int(r["discrete_scale"])
            out_rows.append(
                {
                    "sample": sample,
                    "binary": int(r["binary"]),
                    "discrete_scale": d,
                    "strong_infect": 1 if d == 3 else 0,
                    "weak_infect": 1 if d in {1, 2} else 0,
                    "infection_class": 0 if d == 0 else (1 if d in {1, 2} else 2),
                    "group": infer_group(sample),
                }
            )

    fields = list(out_rows[0].keys()) if out_rows else [
        "sample",
        "binary",
        "discrete_scale",
        "strong_infect",
        "weak_infect",
        "infection_class",
        "group",
    ]
    with OUT_PHENO.open("w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    print(f"Wrote {OUT_PHENO} ({len(out_rows)} rows)")


if __name__ == "__main__":
    main()
