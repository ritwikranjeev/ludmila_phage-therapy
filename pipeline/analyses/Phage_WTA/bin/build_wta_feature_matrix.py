#!/usr/bin/env python3
"""
Build host × WTA feature matrix from raw hmmsearch summaries.

Usage:
    python bin/build_wta_feature_matrix.py

Input:  data/wta_features/raw_wta_hits.tsv (columns: host_id, gene, best_evalue, bitscore, hit_length)
Output: data/wta_features/wta_feature_matrix.tsv
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

# ---------------------------------------------------------------------------
# Presence/absence threshold (full-sequence E-value). Keep in sync with
# WTA_EVALUE_MAX in bin/scan_hosts_for_wta_genes.py.
# ---------------------------------------------------------------------------
WTA_EVALUE_MAX: float = 1e-5

# Gene columns in matrix order (basename without _ref suffix).
GENE_COLUMNS: list[str] = [
    "tarM",
    "tarS",
    "tarP",
    "tagO",
    "tarI",
    "tarJ",
    "tarK",
    "tarL",
    "tagN",
    "tagE",
]


def project_root() -> Path:
    return Path(__file__).resolve().parent.parent


def read_host_ids(tsv_path: Path) -> list[str]:
    if not tsv_path.is_file():
        return []
    ids: list[str] = []
    for line in tsv_path.read_text(encoding="utf-8").splitlines()[1:]:
        line = line.strip()
        if not line:
            continue
        ids.append(line.split("\t", 1)[0].strip())
    return ids


def gene_to_column(gene: str) -> str | None:
    g = gene.strip()
    if g.endswith("_ref"):
        g = g[: -len("_ref")]
    return g if g in GENE_COLUMNS else None


def infer_backbone(present: dict[str, int]) -> str:
    # -----------------------------------------------------------------------
    # TODO: Replace with literature-based rules linking WTA biosynthesis gene
    # profiles (tar/tag loci) to ribitol-phosphate (RboP) vs glycerol-phosphate
    # (GroP) class WTAs and to species-specific variants.
    # -----------------------------------------------------------------------
    if present.get("tarL", 0):
        return "RboP"
    if present.get("tarM", 0) and present.get("tagO", 0):
        return "GroP"
    return "unknown"


def infer_glycan_flags(present: dict[str, int]) -> dict[str, int]:
    # -----------------------------------------------------------------------
    # TODO: Map gene presence (e.g. glycosyltransferases) to predicted
    # alpha-GlcNAc, beta-GlcNAc, GalNAc, Glc decoration flags using curated
    # genotype–phenotype rules from WTA reviews / key papers.
    # -----------------------------------------------------------------------
    return {
        "alpha_GlcNAc": 0,
        "beta_GlcNAc": 0,
        "GalNAc": 0,
        "Glc": 0,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="Build WTA feature matrix from hmmsearch summaries.")
    ap.add_argument(
        "--allow-empty-raw",
        action="store_true",
        help="If raw_wta_hits.tsv is missing, write a zeroed matrix from config/host_genomes.tsv only.",
    )
    args = ap.parse_args()

    root = project_root()
    raw_path = root / "data" / "wta_features" / "raw_wta_hits.tsv"
    out_path = root / "data" / "wta_features" / "wta_feature_matrix.tsv"

    ev: dict[str, dict[str, float]] = defaultdict(dict)
    if not raw_path.is_file():
        if not args.allow_empty_raw:
            print(f"ERROR: Missing {raw_path}", file=sys.stderr)
            return 1
        print(
            f"NOTE: {raw_path} missing; building matrix from host list only (all gene flags 0).",
            file=sys.stderr,
        )
    else:
        with raw_path.open(encoding="utf-8") as fh:
            header = fh.readline()
            if not header.strip():
                print("ERROR: Empty raw hits file.", file=sys.stderr)
                return 1
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                host_id, gene, ev_s, _bitscore, _hitlen = (
                    parts[0],
                    parts[1],
                    parts[2],
                    parts[3],
                    parts[4],
                )
                col = gene_to_column(gene)
                if col is None:
                    continue
                try:
                    evalue = float(ev_s)
                except ValueError:
                    continue
                cur = ev[host_id].get(col)
                if cur is None or evalue < cur:
                    ev[host_id][col] = evalue

    # Cohort = config/host_genomes.tsv only (do not append extra host_ids seen only in raw hits).
    tsv_hosts = read_host_ids(root / "config" / "host_genomes.tsv")
    if tsv_hosts:
        hosts_sorted = list(dict.fromkeys(tsv_hosts))
    else:
        hosts_sorted = sorted(ev.keys())
    gly_keys = ["alpha_GlcNAc", "beta_GlcNAc", "GalNAc", "Glc"]
    out_cols = (
        ["host_id", "backbone_type"]
        + gly_keys
        + GENE_COLUMNS
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(out_cols) + "\n")
        for host_id in hosts_sorted:
            present = {
                g: int(ev[host_id].get(g, float("inf")) <= WTA_EVALUE_MAX)
                for g in GENE_COLUMNS
            }
            backbone = infer_backbone(present)
            gly = infer_glycan_flags(present)
            row = [host_id, backbone]
            row += [str(gly[k]) for k in gly_keys]
            row += [str(present[g]) for g in GENE_COLUMNS]
            fh.write("\t".join(row) + "\n")

    print(f"Wrote {out_path} ({len(hosts_sorted)} host(s)).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
