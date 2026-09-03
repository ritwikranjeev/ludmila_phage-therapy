#!/usr/bin/env python3
"""Create a lenient Medaka-only QC table for manual review."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path
from typing import Dict, Tuple


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


def calc_af(fmt_d: Dict[str, str], info_d: Dict[str, str], dp: float) -> Tuple[str, str]:
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
    dps = info_d.get("DPS", "")
    if dps and "," in dps:
        try:
            nums = [float(x) for x in dps.split(",") if x not in ("", ".")]
            if len(nums) >= 2 and sum(nums) > 0:
                return dps, f"{nums[1] / sum(nums):.4f}"
        except ValueError:
            pass
    ao = fmt_d.get("AO", "") or info_d.get("AO", "")
    try:
        if ao and dp > 0:
            return ad, f"{float(ao) / dp:.4f}"
    except ValueError:
        pass
    return ad, ""


def region_flag(effect: str, chrom: str, pos: int, rrna_window_end: int) -> str:
    e = (effect or "").lower()
    if chrom == "contig_2" and pos <= rrna_window_end:
        return "rRNA_adjacent"
    if "upstream" in e or "downstream" in e:
        return "upstream_downstream"
    if "intergenic" in e:
        return "intergenic"
    return "coding_or_structural"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--min-qual", type=float, default=20.0)  # lenient default
    p.add_argument("--min-dp", type=float, default=10.0)    # lenient default
    p.add_argument("--min-af", type=float, default=0.8)     # lenient default
    p.add_argument("--rrna-window-end", type=int, default=11000)
    p.add_argument("--annotation-tsv", type=Path, default=Path("tables/mutations.tsv"))
    args = p.parse_args()

    ann_map: Dict[Tuple[str, int, str, str], Dict[str, str]] = {}
    if args.annotation_tsv.exists():
        with args.annotation_tsv.open(newline="") as af:
            rr = csv.DictReader(af, delimiter="\t")
            for row in rr:
                try:
                    ann_map[(row["chrom"], int(row["pos"]), row["ref"], row["alt"])] = row
                except Exception:
                    continue

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "qual",
        "filter",
        "gt",
        "dp",
        "ad",
        "af",
        "gene",
        "product",
        "effect",
        "impact",
        "region",
    ]
    with open_text(args.vcf) as fh, args.output.open("w", newline="") as oh:
        w = csv.DictWriter(oh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 10:
                continue
            chrom, pos_s, _id, ref, alt, qual_s, filt, info, fmt, sample = c[:10]
            pos = int(pos_s)
            try:
                qual = float(qual_s) if qual_s not in (".", "") else 0.0
            except ValueError:
                qual = 0.0
            info_d = parse_info(info)
            fmt_d = parse_fmt(fmt, sample)
            gt = fmt_d.get("GT", "")
            try:
                dp = float(fmt_d.get("DP", info_d.get("DP", "0")) or 0)
            except ValueError:
                dp = 0.0
            ad, af_s = calc_af(fmt_d, info_d, dp)
            try:
                af = float(af_s) if af_s else -1.0
            except ValueError:
                af = -1.0
            ann = info_d.get("ANN", "")
            top = ann.split(",")[0] if ann else ""
            bits = top.split("|") if top else []
            effect = bits[1] if len(bits) > 1 else ""
            impact = bits[2] if len(bits) > 2 else ""
            gene = bits[3] if len(bits) > 3 else ""
            product = bits[10] if len(bits) > 10 else ""
            ann_row = ann_map.get((chrom, pos, ref, alt))
            if ann_row:
                gene = gene or ann_row.get("gene", "")
                product = product or ann_row.get("product", "")
                effect = effect or ann_row.get("effect_detail", "") or ann_row.get("effect", "")
                impact = impact or ann_row.get("effect", "")
            q_ok = qual >= args.min_qual
            d_ok = dp >= args.min_dp
            af_ok = True if af < 0 else (af >= args.min_af)
            if not (q_ok and d_ok and af_ok):
                continue
            w.writerow(
                {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "qual": f"{qual:.2f}",
                    "filter": filt,
                    "gt": gt,
                    "dp": f"{dp:.2f}",
                    "ad": ad,
                    "af": af_s,
                    "gene": gene,
                    "product": product,
                    "effect": effect,
                    "impact": impact,
                    "region": region_flag(effect, chrom, pos, args.rrna_window_end),
                }
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
