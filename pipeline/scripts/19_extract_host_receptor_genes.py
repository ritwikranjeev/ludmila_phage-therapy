#!/usr/bin/env python3
"""Scan Prokka GFFs for WTA / surface receptor candidate genes (curated panel).

Writes:
  outputs/receptor_candidates/host_receptor_presence.Rtab — pyseer --pres (Gene × samples, 0/1)
  outputs/receptor_candidates/host_receptor_long.tsv — one row per gene hit per genome
  outputs/receptor_candidates/host_receptor_allele_summary.tsv — copy count + bp span per (sample, gene)

Allele detail is coarse (n copies, bp lengths) without variant calling; extend with VCF later if needed.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PANEL_PATH = PROJECT_ROOT / "config" / "host_receptor_gene_panel.tsv"
ANNOT_ROOT = PROJECT_ROOT / "outputs" / "annotation"
OUT_DIR = PROJECT_ROOT / "outputs" / "receptor_candidates"
SAMPLES_TSV = PROJECT_ROOT / "config" / "samples.tsv"


def load_panel(path: Path) -> Tuple[Dict[str, str], Dict[str, List[str]]]:
    """Return canonical_id -> tier, and synonym_lower -> canonical_id."""
    tier: Dict[str, str] = {}
    synmap: Dict[str, str] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            cid = row["canonical_id"].strip()
            tier[cid] = row["tier"].strip()
            syns = [s.strip() for s in row["synonyms"].split(",") if s.strip()]
            for s in syns:
                synmap[s.lower()] = cid
            synmap[cid.lower()] = cid
    return tier, synmap


def parse_gff_attrs(field9: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for part in field9.split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        k, v = part.split("=", 1)
        out[k.strip().lower()] = v.strip()
    return out


def iter_cds_features(gff_path: Path) -> Iterator[Dict[str, object]]:
    with gff_path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            if cols[2] != "CDS":
                continue
            attrs = parse_gff_attrs(cols[8])
            gene = (attrs.get("gene") or "").strip()
            yield {
                "contig": cols[0],
                "start": int(cols[3]),
                "end": int(cols[4]),
                "strand": cols[6],
                "gene": gene,
                "locus_tag": (attrs.get("locus_tag") or "").strip(),
                "product": (attrs.get("product") or "").strip(),
            }


def match_canonical(gene_field: str, synmap: Dict[str, str]) -> Optional[str]:
    if not gene_field:
        return None
    return synmap.get(gene_field.lower())


def tier_rank(tier: str) -> int:
    order = [
        "core_WTA",
        "WTA_glyco",
        "WTA_tailoring",
        "regulatory_two_component",
        "regulatory_cell_wall_stress",
        "regulatory_quorum_sensing",
    ]
    try:
        return order.index(tier)
    except ValueError:
        return 99


def main() -> None:
    tier_map, synmap = load_panel(PANEL_PATH)
    samples: List[str] = []
    with SAMPLES_TSV.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            samples.append(row["sample"].strip())

    presence: Dict[str, Dict[str, int]] = defaultdict(lambda: {s: 0 for s in samples})
    long_rows: List[dict] = []
    allele_rows: List[dict] = []

    for sample in samples:
        gff = ANNOT_ROOT / sample / f"{sample}.gff"
        hits: Dict[str, List[dict]] = defaultdict(list)
        if gff.is_file():
            for feat in iter_cds_features(gff):
                cid = match_canonical(str(feat["gene"]), synmap)
                if cid:
                    hits[cid].append(feat)
        for cid in tier_map:
            feats = hits.get(cid, [])
            presence[cid][sample] = 1 if feats else 0
            span_bp = sum(abs(f["end"] - f["start"]) + 1 for f in feats) if feats else 0
            allele_rows.append(
                {
                    "sample": sample,
                    "canonical_id": cid,
                    "tier": tier_map[cid],
                    "n_copies": len(feats),
                    "span_bp_sum": span_bp,
                    "locus_tags": ",".join(str(f.get("locus_tag") or "") for f in feats),
                    "contigs": ",".join(str(f["contig"]) for f in feats),
                }
            )
            for f in feats:
                long_rows.append(
                    {
                        "sample": sample,
                        "canonical_id": cid,
                        "tier": tier_map[cid],
                        "locus_tag": f.get("locus_tag"),
                        "contig": f["contig"],
                        "start": f["start"],
                        "end": f["end"],
                        "strand": f["strand"],
                        "product": f.get("product"),
                    }
                )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    gene_ids = sorted(
        tier_map.keys(),
        key=lambda g: (tier_rank(tier_map[g]), g),
    )

    rtab_path = OUT_DIR / "host_receptor_presence.Rtab"
    with rtab_path.open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["Gene"] + samples)
        for gid in gene_ids:
            row = [f"host_{gid}"] + [str(presence[gid][s]) for s in samples]
            w.writerow(row)
    print(f"Wrote {rtab_path}")

    long_path = OUT_DIR / "host_receptor_long.tsv"
    fields = [
        "sample",
        "canonical_id",
        "tier",
        "locus_tag",
        "contig",
        "start",
        "end",
        "strand",
        "product",
    ]
    with long_path.open("w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(long_rows)
    print(f"Wrote {long_path} ({len(long_rows)} hit rows)")

    sum_path = OUT_DIR / "host_receptor_allele_summary.tsv"
    af = [
        "sample",
        "canonical_id",
        "tier",
        "n_copies",
        "span_bp_sum",
        "locus_tags",
        "contigs",
    ]
    with sum_path.open("w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=af, delimiter="\t")
        w.writeheader()
        w.writerows(allele_rows)
    print(f"Wrote {sum_path}")


if __name__ == "__main__":
    main()
