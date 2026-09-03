#!/usr/bin/env python3
"""Combine SnpEff-annotated VCFs into one TSV mutation table."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("r")


def parse_vcf_samples(header_line: str) -> List[str]:
    if not header_line.startswith("#CHROM"):
        return []
    return header_line.strip().split("\t")[9:]


def split_info_ann(ann_raw: str) -> List[str]:
    anns = []
    for part in ann_raw.split(","):
        part = part.strip()
        if part.startswith("ANN="):
            part = part[4:]
        fields = part.split("|")
        if len(fields) >= 5:
            anns.append(part)
    return anns


def pick_main_ann(ann_list: List[str]) -> str:
    """Prefer a usable SnpEff annotation (skip ERROR_CHROMOSOME_NOT_FOUND stubs)."""
    for ann in ann_list:
        if ann and "ERROR_CHROMOSOME" not in ann:
            return ann
    return ann_list[0] if ann_list else ""


def classify_effect(annotation: str) -> str:
    a = annotation.lower()
    if "frameshift" in a:
        return "frameshift"
    if "stop_gained" in a or "stop gained" in a:
        return "stop-gain"
    if "start_lost" in a:
        return "start-lost"
    if "missense" in a:
        return "nonsynonymous"
    if "synonymous" in a:
        return "synonymous"
    if "intergenic" in a:
        return "intergenic"
    if "upstream" in a or "downstream" in a:
        return "regulatory"
    return "other"


def parse_ann_fields(ann: str) -> Tuple[str, str, str, str]:
    bits = ann.split("|")
    effect = bits[1] if len(bits) > 1 else ""
    gene = bits[3] if len(bits) > 3 else ""
    gene_id = bits[4] if len(bits) > 4 else ""
    product = bits[10] if len(bits) > 10 else ""
    return effect, gene or gene_id, product, ann


def iter_vcf_records(path: Path) -> Iterable[Tuple[dict, str]]:
    with open_text(path) as fh:
        header_samples: List[str] = []
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header_samples = parse_vcf_samples(line)
                continue
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, _id, ref, alt, qual, filt, info = cols[:8]
            fmt = cols[8] if len(cols) > 8 else ""
            samples = cols[9:] if len(cols) > 9 else []
            yield {
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "qual": qual,
                "filter": filt,
                "info": info,
                "format": fmt,
                "samples": samples,
                "header_samples": header_samples,
            }, line


def variant_key(rec: dict) -> Tuple[str, str, str, str]:
    return (rec["chrom"], rec["pos"], rec["ref"], rec["alt"])


def sniff_variant_type(ref: str, alt: str) -> str:
    alts = [a for a in alt.split(",") if a]
    if len(alts) != 1:
        return "multi"
    a = alts[0]
    if a.startswith("<") and a.endswith(">"):
        return "SV"
    if len(ref) == 1 and len(a) == 1:
        return "SNP"
    return "indel"


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--vcf",
        action="append",
        required=True,
        metavar="NAME:PATH",
        help="Sample label and VCF path (repeatable), e.g. pre:vcf/pre.snv.snpeff.vcf.gz",
    )
    p.add_argument("-o", "--output", required=True, help="Output TSV path")
    p.add_argument(
        "--reference",
        required=True,
        metavar="NAME",
        help="Sample ID treated as reference (e.g. pre-phage assembly / ref genome)",
    )
    p.add_argument(
        "--query",
        required=True,
        metavar="NAME",
        help="Sample ID compared to reference (e.g. post-phage)",
    )
    args = p.parse_args(argv)

    labeled: List[Tuple[str, Path]] = []
    for spec in args.vcf:
        if ":" not in spec:
            print("Expected --vcf sample:path", file=sys.stderr)
            return 2
        name, path = spec.split(":", 1)
        labeled.append((name, Path(path)))

    seen_labels = {n for n, _ in labeled}
    if args.reference not in seen_labels:
        print(f"Reference sample {args.reference!r} not in --vcf list", file=sys.stderr)
        return 2
    if args.query not in seen_labels:
        print(f"Query sample {args.query!r} not in --vcf list", file=sys.stderr)
        return 2

    by_key: Dict[Tuple[str, str, str, str], Dict[str, dict]] = defaultdict(dict)

    for sample_name, vcf_path in labeled:
        if not vcf_path.is_file():
            print(f"Missing VCF: {vcf_path}", file=sys.stderr)
            return 2
        for rec, _raw in iter_vcf_records(vcf_path):
            key = variant_key(rec)
            info = rec["info"]
            ann_match = re.search(r"ANN=([^;\t]+)", info)
            ann_raw = ann_match.group(1) if ann_match else ""
            ann_list = split_info_ann(ann_raw)
            main_ann = pick_main_ann(ann_list)
            eff, gene, product, ann_full = parse_ann_fields(main_ann) if main_ann else ("", "", "", "")
            vtype = sniff_variant_type(rec["ref"], rec["alt"])
            effect_bucket = classify_effect(main_ann) if main_ann else "unknown"

            by_key[key][sample_name] = {
                "chrom": rec["chrom"],
                "pos": rec["pos"],
                "ref": rec["ref"],
                "alt": rec["alt"],
                "variant_type": vtype,
                "effect_detail": eff,
                "effect": effect_bucket,
                "gene": gene,
                "product": product,
                "annotation": ann_full,
                "qual": rec["qual"],
                "filter": rec["filter"],
            }

    ref = args.reference
    qry = args.query
    col_in_ref = f"present_in_{ref}"
    col_in_qry = f"present_in_{qry}"
    col_unique_qry = f"unique_to_{qry}"

    base_fields = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "variant_type",
        "effect",
        "effect_detail",
        "gene",
        "product",
        "annotation",
        "qual",
        "filter",
    ]
    fields = base_fields + [col_in_ref, col_in_qry, col_unique_qry]

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    def sort_key(k: Tuple[str, str, str, str]):
        chrom, pos, _r, _a = k
        try:
            ip = int(pos)
        except ValueError:
            ip = 0
        return chrom, ip

    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for key in sorted(by_key.keys(), key=sort_key):
            rows = by_key[key]
            base = next(iter(rows.values()))
            in_ref = ref in rows
            in_qry = qry in rows
            unique_qry = in_qry and not in_ref
            row = {f: base.get(f, "") for f in base_fields}
            row[col_in_ref] = str(in_ref)
            row[col_in_qry] = str(in_qry)
            row[col_unique_qry] = str(unique_qry)
            w.writerow(row)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
