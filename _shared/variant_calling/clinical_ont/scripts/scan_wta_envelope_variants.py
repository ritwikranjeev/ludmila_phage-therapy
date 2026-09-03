#!/usr/bin/env python3
"""Scan high-confidence variants for WTA and envelope-related hits."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import yaml


def load_catalog(path: Path):
    data = yaml.safe_load(path.read_text()) or {}
    categories = data.get("categories", {}) or {}
    product_keywords = data.get("product_keywords", {}) or {}
    gene_to_category = {}
    for cat, genes in categories.items():
        for g in genes or []:
            gene_to_category[str(g).lower()] = str(cat)
    return gene_to_category, product_keywords


def load_focus_tags(path: Path) -> set[str]:
    out = set()
    if not path.exists():
        return out
    with path.open() as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            out.add(s)
    return out


def pick_category(gene: str, product: str, gene_to_cat, product_keywords) -> str:
    g = (gene or "").lower()
    p = (product or "").lower()
    if g in gene_to_cat:
        return gene_to_cat[g]
    for cat, words in product_keywords.items():
        for w in words or []:
            if str(w).lower() in p:
                return str(cat)
    return "other"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--summary", type=Path, required=True)
    p.add_argument("--catalog", type=Path, required=True)
    p.add_argument("--focus-locus-tags", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()

    gene_to_cat, product_keywords = load_catalog(args.catalog)
    focus_tags = load_focus_tags(args.focus_locus_tags)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.summary.open(newline="") as fh, args.output.open("w", newline="") as oh:
        r = csv.DictReader(fh, delimiter="\t")
        fields = [
            "chrom",
            "pos",
            "ref",
            "alt",
            "locus_tag",
            "gene",
            "product",
            "effect",
            "impact",
            "category",
            "tier",
            "caller_support",
            "pass_both_callers",
            "region_flag",
        ]
        w = csv.DictWriter(oh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for row in r:
            locus_tag = row.get("locus_tag", "")
            gene = row.get("gene", "")
            product = row.get("product", "")
            category = pick_category(gene, product, gene_to_cat, product_keywords)
            keep = (category != "other") or (locus_tag in focus_tags)
            if not keep:
                continue
            out = {k: row.get(k, "") for k in fields}
            out["category"] = category if category != "other" else row.get("category", "other")
            w.writerow(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
