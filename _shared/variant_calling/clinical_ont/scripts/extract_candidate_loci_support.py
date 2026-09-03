#!/usr/bin/env python3
"""Extract per-caller support for candidate loci from annotated VCFs."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("r")


def parse_info(info: str) -> Dict[str, str]:
    d: Dict[str, str] = {}
    for p in info.split(";"):
        if not p:
            continue
        if "=" in p:
            k, v = p.split("=", 1)
            d[k] = v
        else:
            d[p] = "True"
    return d


def parse_fmt(fmt: str, sample: str) -> Dict[str, str]:
    k = fmt.split(":") if fmt else []
    v = sample.split(":") if sample else []
    return {kk: (v[i] if i < len(v) else "") for i, kk in enumerate(k)}


def parse_ad_af(fmt_d: Dict[str, str], info_d: Dict[str, str], dp: float) -> Tuple[str, str]:
    ad = fmt_d.get("AD", "") or info_d.get("AD", "")
    af = fmt_d.get("AF", "") or info_d.get("AF", "")
    if af:
        return ad, af
    if ad and "," in ad:
        try:
            nums = [float(x) for x in ad.split(",") if x not in ("", ".")]
            if len(nums) >= 2 and sum(nums) > 0:
                return ad, f"{nums[1] / sum(nums):.4f}"
        except ValueError:
            pass
    if dp > 0:
        ao = fmt_d.get("AO", "") or info_d.get("AO", "")
        try:
            if ao:
                return ad, f"{float(ao) / dp:.4f}"
        except ValueError:
            pass
    return ad, ""


def candidate_hit(gene: str, product: str, chrom: str, pos: int, terms: List[str], win: int = 500) -> bool:
    g = (gene or "").lower()
    p = (product or "").lower()
    for t in terms:
        tl = t.lower()
        if tl in g or tl in p:
            return True
    # explicit cdaR upstream window on contig_2 commonly discussed
    if "cdar_window" in [x.lower() for x in terms] and chrom == "contig_2":
        # user can refine this if exact cdaR coordinate differs
        return False
    return False


def load_annotation_table(path: Path) -> Dict[Tuple[str, int, str, str], Dict[str, str]]:
    out: Dict[Tuple[str, int, str, str], Dict[str, str]] = {}
    if not path.exists():
        return out
    with path.open(newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            try:
                out[(row["chrom"], int(row["pos"]), row["ref"], row["alt"])] = row
            except Exception:
                continue
    return out


def iter_rows(vcf: Path, caller: str, terms: List[str], ann_table: Dict[Tuple[str, int, str, str], Dict[str, str]]) -> Iterable[Dict[str, str]]:
    with open_text(vcf) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 10:
                continue
            chrom, pos, _id, ref, alt, qual, filt, info, fmt, sample = c[:10]
            info_d = parse_info(info)
            ann = info_d.get("ANN", "")
            top = ann.split(",")[0] if ann else ""
            bits = top.split("|") if top else []
            effect = bits[1] if len(bits) > 1 else ""
            gene = bits[3] if len(bits) > 3 else ""
            product = bits[10] if len(bits) > 10 else ""
            ann_row = ann_table.get((chrom, int(pos), ref, alt))
            if ann_row:
                gene = gene or ann_row.get("gene", "")
                product = product or ann_row.get("product", "")
                effect = effect or ann_row.get("effect_detail", "") or ann_row.get("effect", "")
            if not candidate_hit(gene, product, chrom, int(pos), terms):
                continue
            fmt_d = parse_fmt(fmt, sample)
            gt = fmt_d.get("GT", "")
            try:
                dp = float(fmt_d.get("DP", info_d.get("DP", "0")) or 0)
            except ValueError:
                dp = 0.0
            ad, af = parse_ad_af(fmt_d, info_d, dp)
            yield {
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "caller": caller,
                "qual": qual,
                "filter": filt,
                "gt": gt,
                "dp": f"{dp:.2f}",
                "ad": ad,
                "af": af,
                "gene": gene,
                "product": product,
                "effect": effect,
            }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--medaka-vcf", type=Path, required=True)
    p.add_argument("--clair3-vcf", type=Path, required=True)
    p.add_argument("--annotation-tsv", type=Path, default=Path("tables/mutations.tsv"))
    p.add_argument("--terms", default="agrA,cdaR,mdrP,ebh,rpsC")
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()

    terms = [x.strip() for x in args.terms.split(",") if x.strip()]
    ann_table = load_annotation_table(args.annotation_tsv)
    rows = list(iter_rows(args.medaka_vcf, "Medaka", terms, ann_table)) + list(
        iter_rows(args.clair3_vcf, "Clair3", terms, ann_table)
    )
    rows.sort(key=lambda r: (r["chrom"], int(r["pos"]), r["caller"]))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "caller",
        "qual",
        "filter",
        "gt",
        "dp",
        "ad",
        "af",
        "gene",
        "product",
        "effect",
    ]
    with args.output.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
