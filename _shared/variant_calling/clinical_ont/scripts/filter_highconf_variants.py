#!/usr/bin/env python3
"""Build a high-confidence Medaka+Clair3 variant set with simple QC metadata."""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import yaml


@dataclass
class VcfRec:
    chrom: str
    pos: int
    ref: str
    alt: str
    qual: float
    filt: str
    info: str
    fmt: str
    sample: str
    raw: str


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("r")


def parse_info(info: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for part in info.split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
        else:
            out[part] = "True"
    return out


def parse_gt_sample(fmt: str, sample: str) -> Dict[str, str]:
    keys = fmt.split(":") if fmt else []
    vals = sample.split(":") if sample else []
    out: Dict[str, str] = {}
    for i, k in enumerate(keys):
        out[k] = vals[i] if i < len(vals) else ""
    return out


def first_ann(info: str) -> Tuple[str, str, str, str]:
    fields = parse_info(info)
    ann = fields.get("ANN", "")
    if not ann:
        return "", "", "", ""
    top = ann.split(",")[0]
    bits = top.split("|")
    effect = bits[1] if len(bits) > 1 else ""
    impact = bits[2] if len(bits) > 2 else ""
    gene = bits[3] if len(bits) > 3 else ""
    product = bits[10] if len(bits) > 10 else ""
    return effect, impact, gene, product


def classify_effect(effect: str) -> str:
    e = effect.lower()
    if "missense" in e:
        return "missense"
    if "synonymous" in e:
        return "synonymous"
    if "upstream" in e:
        return "upstream"
    if "downstream" in e:
        return "downstream"
    if "intergenic" in e:
        return "intergenic"
    if "frameshift" in e:
        return "frameshift"
    if "stop_gained" in e:
        return "stop_gained"
    if "intron" in e:
        return "intron"
    return effect or "other"


def parse_vcf_records(path: Path) -> Tuple[List[str], Dict[Tuple[str, int, str, str], VcfRec]]:
    headers: List[str] = []
    recs: Dict[Tuple[str, int, str, str], VcfRec] = {}
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##"):
                headers.append(line.rstrip("\n"))
                continue
            if line.startswith("#CHROM"):
                headers.append(line.rstrip("\n"))
                continue
            if not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 10:
                continue
            qual = 0.0
            try:
                qual = float(c[5]) if c[5] not in (".", "") else 0.0
            except ValueError:
                qual = 0.0
            rec = VcfRec(
                chrom=c[0],
                pos=int(c[1]),
                ref=c[3],
                alt=c[4],
                qual=qual,
                filt=c[6],
                info=c[7],
                fmt=c[8],
                sample=c[9],
                raw=line.rstrip("\n"),
            )
            recs[(rec.chrom, rec.pos, rec.ref, rec.alt)] = rec
    return headers, recs


def parse_vcf_records_by_site(path: Path) -> Dict[Tuple[str, int, str], List[VcfRec]]:
    out: Dict[Tuple[str, int, str], List[VcfRec]] = {}
    _, recs = parse_vcf_records(path)
    for rec in recs.values():
        out.setdefault((rec.chrom, rec.pos, rec.ref), []).append(rec)
    return out


def parse_bed(path: Path) -> Dict[str, List[Tuple[int, int, str]]]:
    out: Dict[str, List[Tuple[int, int, str]]] = {}
    if not path.exists():
        return out
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            label = cols[3] if len(cols) > 3 else "hard_region"
            out.setdefault(chrom, []).append((start, end, label))
    return out


def region_flag(chrom: str, pos1: int, regions: Dict[str, List[Tuple[int, int, str]]]) -> str:
    for start0, end1, label in regions.get(chrom, []):
        if start0 <= (pos1 - 1) < end1:
            return label
    return "other"


def parse_catalog(path: Path) -> Dict[str, str]:
    data = yaml.safe_load(path.read_text()) or {}
    out: Dict[str, str] = {}
    for cat, genes in (data.get("categories") or {}).items():
        for g in genes or []:
            out[str(g).lower()] = str(cat)
    return out


def load_focus_tags(path: Path) -> set[str]:
    tags = set()
    if not path.exists():
        return tags
    with path.open() as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            tags.add(s)
    return tags


def pick_num(fmt_map: Dict[str, str], keys: List[str], default: float = 0.0) -> float:
    for k in keys:
        v = fmt_map.get(k, "")
        if not v or v == ".":
            continue
        # AD like ref,alt
        if "," in v and k == "AD":
            try:
                a = [float(x) for x in v.split(",") if x not in ("", ".")]
                if len(a) >= 2 and sum(a) > 0:
                    return a[1] / sum(a)
            except ValueError:
                pass
            continue
        try:
            return float(v)
        except ValueError:
            continue
    return default


def get_dp_af(rec: VcfRec) -> Tuple[float, float]:
    info = parse_info(rec.info)
    fmt_map = parse_gt_sample(rec.fmt, rec.sample)
    dp = pick_num(fmt_map, ["DP"], pick_num(info, ["DP"], 0.0))
    # Prefer explicit AF, else derive from AD, else from AO/DP where present.
    af = pick_num(fmt_map, ["AF"], pick_num(info, ["AF"], -1.0))
    if af < 0:
        af = pick_num(fmt_map, ["AD"], -1.0)
    if af < 0:
        dps = info.get("DPS", "")
        if dps and "," in dps:
            try:
                vals = [float(x) for x in dps.split(",") if x not in ("", ".")]
                if len(vals) >= 2 and sum(vals) > 0:
                    af = vals[1] / sum(vals)
            except ValueError:
                pass
    if af < 0:
        ao = pick_num(fmt_map, ["AO"], pick_num(info, ["AO"], -1.0))
        if ao >= 0 and dp > 0:
            af = ao / dp
    if af < 0:
        af = 0.0
    return dp, af


def load_annotation_table(path: Path) -> Dict[Tuple[str, int, str, str], Dict[str, str]]:
    out: Dict[Tuple[str, int, str, str], Dict[str, str]] = {}
    if not path.exists():
        return out
    with path.open(newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            try:
                key = (row["chrom"], int(row["pos"]), row["ref"], row["alt"])
            except Exception:
                continue
            out[key] = row
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--medaka-vcf", type=Path, required=True)
    p.add_argument("--clair3-vcf", type=Path, required=True)
    p.add_argument("--hard-regions", type=Path, required=True)
    p.add_argument("--focus-catalog", type=Path, required=True)
    p.add_argument("--focus-locus-tags", type=Path, required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--annotation-tsv", type=Path, default=Path("tables/mutations.tsv"))
    # Defaults chosen to include validated controls: agrA, fumC-upstream, yclM.
    p.add_argument("--medaka-min-qual", type=float, default=3.0)
    p.add_argument("--medaka-min-dp", type=float, default=10.0)
    p.add_argument("--clair3-min-qual", type=float, default=8.0)
    p.add_argument("--clair3-min-dp", type=float, default=10.0)
    p.add_argument("--medaka-solo-min-qual", type=float, default=30.0)
    p.add_argument("--medaka-solo-min-dp", type=float, default=15.0)
    p.add_argument("--medaka-solo-min-af", type=float, default=0.9)
    p.add_argument("--allow-medaka-only-in-hard-regions", default="true")
    p.add_argument(
        "--positive-controls",
        default="contig_2:848821:G:C,contig_2:1486788:TA:T,contig_2:1646692:G:A",
    )
    p.add_argument("--output-vcf", type=Path, required=True)
    p.add_argument("--output-summary", type=Path, required=True)
    args = p.parse_args()

    allow_medaka_only = str(args.allow_medaka_only_in_hard_regions).lower() in {
        "true",
        "1",
        "yes",
        "y",
    }
    _, medaka = parse_vcf_records(args.medaka_vcf)
    clair3_by_site = parse_vcf_records_by_site(args.clair3_vcf)
    regions = parse_bed(args.hard_regions)
    gene_to_cat = parse_catalog(args.focus_catalog)
    focus_tags = load_focus_tags(args.focus_locus_tags)
    ann_table = load_annotation_table(args.annotation_tsv)

    args.output_vcf.parent.mkdir(parents=True, exist_ok=True)
    args.output_summary.parent.mkdir(parents=True, exist_ok=True)

    header = [
        "##fileformat=VCFv4.2",
        "##INFO=<ID=PASS_BOTH_CALLERS,Number=1,Type=String,Description=\"true if Medaka and Clair3 match\">",
        "##INFO=<ID=CALLER_SUPPORT,Number=1,Type=String,Description=\"both or medaka_only\">",
        "##INFO=<ID=REGION_FLAG,Number=1,Type=String,Description=\"hard region label or other\">",
        "##INFO=<ID=TIER,Number=1,Type=String,Description=\"Tier 1/2/3 interpretation tag\">",
        "##INFO=<ID=CATEGORY,Number=1,Type=String,Description=\"WTA/envelope category\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ]

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
        "medaka_gt",
        "medaka_dp",
        "medaka_af",
        "medaka_qual",
        "clair3_gt",
        "clair3_dp",
        "clair3_af",
        "clair3_qual",
        "passes_filters",
    ]

    passed_keys: set[Tuple[str, int, str, str]] = set()
    with args.output_vcf.open("w") as vo, args.output_summary.open("w", newline="") as to:
        for h in header:
            vo.write(h + "\n")

        w = csv.DictWriter(to, fieldnames=fields, delimiter="\t")
        w.writeheader()

        for key in sorted(medaka.keys(), key=lambda k: (k[0], k[1], k[2], k[3])):
            m = medaka[key]
            cands = clair3_by_site.get((m.chrom, m.pos, m.ref), [])
            # Primary concordance: same allele at same position (GT may differ in haploid-like bacterial runs).
            c = next((x for x in cands if x.alt == m.alt), None)

            mfmt = parse_gt_sample(m.fmt, m.sample)
            m_gt = mfmt.get("GT", "")
            m_dp, m_af = get_dp_af(m)
            m_qual = m.qual

            if c:
                cfmt = parse_gt_sample(c.fmt, c.sample)
                c_gt = cfmt.get("GT", "")
                c_dp, c_af = get_dp_af(c)
                c_qual = c.qual
            else:
                c_gt, c_dp, c_af, c_qual = "", 0.0, 0.0, 0.0

            pass_both = bool(c) and (m.alt == c.alt)
            rflag = region_flag(m.chrom, m.pos, regions)
            caller_support = "both" if pass_both else "medaka_only"

            q_ok = m_qual >= args.medaka_min_qual
            dp_ok = m_dp >= args.medaka_min_dp
            base_ok = q_ok and dp_ok and m.filt in ("PASS", ".")
            if c:
                clair3_ok = (c_qual >= args.clair3_min_qual) and (c_dp >= args.clair3_min_dp)
                pass_logic = pass_both and clair3_ok
            else:
                solo_af_ok = True if m_af <= 0 else (m_af >= args.medaka_solo_min_af)
                solo_ok = (m_qual >= args.medaka_solo_min_qual) and (m_dp >= args.medaka_solo_min_dp) and solo_af_ok
                # Medaka-only calls allowed only in non-hard regions.
                pass_logic = allow_medaka_only and (rflag == "other") and solo_ok
            passes = base_ok and pass_logic
            if not passes:
                continue
            passed_keys.add((m.chrom, m.pos, m.ref, m.alt))

            effect, impact, gene, product = first_ann(m.info)
            ann_row = ann_table.get((m.chrom, m.pos, m.ref, m.alt))
            if ann_row:
                gene = gene or ann_row.get("gene", "")
                product = product or ann_row.get("product", "")
                effect = effect or ann_row.get("effect_detail", "") or ann_row.get("effect", "")
                impact = impact or ann_row.get("effect", "")
            locus_tag = gene if gene.startswith("DEBJFDCL_") else ""
            category = gene_to_cat.get((gene or "").lower(), "other")
            if category == "other":
                for g_lower, cat in gene_to_cat.items():
                    if g_lower and g_lower in (product or "").lower():
                        category = cat
                        break

            tier = "Tier 3"
            key_gene = (gene or "").lower()
            is_focus_gene = key_gene in {"agra", "cdar", "mdrp", "ebh", "rpsc"} or locus_tag in focus_tags
            if is_focus_gene or category != "other":
                tier = "Tier 1"
            elif classify_effect(effect) in {"missense", "frameshift", "stop_gained"}:
                tier = "Tier 2"

            info_extra = (
                f"PASS_BOTH_CALLERS={'true' if pass_both else 'false'};"
                f"CALLER_SUPPORT={caller_support};REGION_FLAG={rflag};"
                f"TIER={tier.replace(' ', '_')};CATEGORY={category}"
            )
            cols = m.raw.split("\t")
            cols[7] = cols[7] + ";" + info_extra if cols[7] else info_extra
            vo.write("\t".join(cols) + "\n")

            w.writerow(
                {
                    "chrom": m.chrom,
                    "pos": m.pos,
                    "ref": m.ref,
                    "alt": m.alt,
                    "locus_tag": locus_tag,
                    "gene": gene,
                    "product": product,
                    "effect": classify_effect(effect),
                    "impact": impact,
                    "category": category,
                    "tier": tier,
                    "caller_support": caller_support,
                    "pass_both_callers": str(pass_both),
                    "region_flag": rflag,
                    "medaka_gt": m_gt,
                    "medaka_dp": f"{m_dp:.2f}",
                    "medaka_af": f"{m_af:.4f}",
                    "medaka_qual": f"{m_qual:.2f}",
                    "clair3_gt": c_gt,
                    "clair3_dp": f"{c_dp:.2f}" if c else "",
                    "clair3_af": f"{c_af:.4f}" if c else "",
                    "clair3_qual": f"{c_qual:.2f}" if c else "",
                    "passes_filters": "True",
                }
            )
    # Positive-control sanity check with raw-record dump to stderr if missing.
    pcs: List[Tuple[str, int, str, str]] = []
    for tok in [x.strip() for x in args.positive_controls.split(",") if x.strip()]:
        try:
            cchrom, cpos, cref, calt = tok.split(":")
            pcs.append((cchrom, int(cpos), cref, calt))
        except ValueError:
            continue
    for pc in pcs:
        if pc in passed_keys:
            continue
        chrom, pos, ref, alt = pc
        print(
            f"[WARN] Positive-control missing from high-confidence: {chrom}:{pos} {ref}>{alt}",
            file=sys.stderr,
        )
        m = medaka.get((chrom, pos, ref, alt))
        if m:
            print(f"[WARN] Medaka raw: {m.raw}", file=sys.stderr)
        else:
            print("[WARN] Medaka raw: <not found>", file=sys.stderr)
        cands = clair3_by_site.get((chrom, pos, ref), [])
        if not cands:
            print("[WARN] Clair3 raw: <not found at site>", file=sys.stderr)
        else:
            for c in cands:
                print(f"[WARN] Clair3 raw: {c.raw}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
