#!/usr/bin/env python3
"""Copy key pipeline outputs into outputs/analysis_planning_bundle/ for offline planning."""

import shutil
from pathlib import Path
from datetime import datetime, timezone

PROJECT_ROOT = Path(__file__).resolve().parents[1]
BUNDLE = PROJECT_ROOT / "outputs" / "analysis_planning_bundle"


def cp(src: Path, dest: Path, manifest: list) -> None:
    if not src.is_file():
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)
    manifest.append((src.relative_to(PROJECT_ROOT), dest.relative_to(BUNDLE)))


def main() -> int:
    manifest = []
    if BUNDLE.exists():
        shutil.rmtree(BUNDLE)
    BUNDLE.mkdir(parents=True)

    # --- GWAS ---
    fdr = PROJECT_ROOT / "outputs" / "gwas" / "post_gwas_fdr"
    if fdr.is_dir():
        for p in sorted(fdr.glob("*.tsv")):
            cp(p, BUNDLE / "gwas" / "fdr_tables" / p.name, manifest)
    cp(
        PROJECT_ROOT / "outputs" / "gwas" / "gwas_priority_overview.tsv",
        BUNDLE / "gwas" / "gwas_priority_overview.tsv",
        manifest,
    )
    for sub in (
        "pyseer_gene",
        "pyseer_defense",
        "pyseer_prophage",
        "pyseer_host_receptors",
    ):
        th = PROJECT_ROOT / "outputs" / "gwas" / sub / "top_hits.tsv"
        cp(th, BUNDLE / "gwas" / "top_hits" / f"{sub}_top_hits.tsv", manifest)
    cp(
        PROJECT_ROOT / "outputs" / "gwas" / "pyseer_gene" / "top_hits_binary_full_fdr.tsv",
        BUNDLE / "gwas" / "top_hits" / "pyseer_gene_top_hits_binary_full_fdr.tsv",
        manifest,
    )

    # --- Integration (phage–host) ---
    integ = PROJECT_ROOT / "outputs" / "integration"
    if integ.is_dir():
        for p in sorted(integ.glob("*.tsv")):
            cp(p, BUNDLE / "integration" / p.name, manifest)

    # --- Phage (parsed) ---
    for name in ("phage_genes.tsv", "phage_receptor_associated.tsv"):
        cp(PROJECT_ROOT / "data" / "phage" / name, BUNDLE / "phage" / name, manifest)

    # --- Config / cohort inputs used in GWAS & integration ---
    for rel in (
        "config/phenotypes.tsv",
        "config/samples.tsv",
        "config/host_receptor_gene_panel.tsv",
        "data/phenotypes/host_range.tsv",
    ):
        src = PROJECT_ROOT / rel
        cp(src, BUNDLE / "config_and_phenotypes" / Path(rel).name, manifest)

    # --- Population structure (if present) ---
    cp(
        PROJECT_ROOT / "outputs" / "gwas" / "population_structure_dist.tsv",
        BUNDLE / "gwas" / "population_structure_dist.tsv",
        manifest,
    )

    # --- Optional summaries (if synced from HPC) ---
    for rel in (
        "outputs/prophages/prophage_summary.tsv",
        "outputs/prophages/prophage_feature_matrix.tsv",
        "outputs/defense/defense_family_presence.tsv",
        "outputs/defense/defense_islands.tsv",
        "outputs/trees/core.treefile",
    ):
        src = PROJECT_ROOT / rel
        sub = Path(rel).parts[1]  # prophages, defense, trees
        cp(src, BUNDLE / sub / Path(rel).name, manifest)

    # Manifest
    lines = [
        f"# Analysis planning bundle",
        f"# Generated: {datetime.now(timezone.utc).isoformat()}",
        f"# Project: {PROJECT_ROOT}",
        "",
        "relative_source -> bundle_path",
        "",
    ]
    for a, b in manifest:
        lines.append(f"{a}\t->\t{b}")
    (BUNDLE / "MANIFEST.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    readme = BUNDLE / "README.txt"
    readme.write_text(
        """This folder is a **snapshot copy** of outputs useful for planning the next analysis stage.

Layout
------
gwas/fdr_tables/     BH-FDR tables (one per analysis x trait) from scripts/run_post_gwas_fdr.sh
gwas/top_hits/       Trimmed tables for notebooks / quick review (binary trait)
gwas/gwas_priority_overview.tsv   Merged preview (top rows per analysis; exploratory only)
gwas/population_structure_dist.tsv   Patristic distance matrix for pyseer (if present under outputs/gwas/)

integration/         Phage–host crosswalk + host receptor GWAS snippet + panel copy

phage/               Parsed Ludmila CDS tables from scripts/17_parse_phage_annotation.py

config_and_phenotypes/   samples.tsv, phenotypes.tsv, host_range.tsv, host_receptor_gene_panel.tsv

prophages/ defense/ trees/   Copied only if present under outputs/ (run this script on HPC to pick them up)

MANIFEST.txt         Exact source -> destination mapping for this run.

Refresh
-------
From project root:
  python3 scripts/build_analysis_planning_bundle.py

If large files are missing, rsync the corresponding paths from HPC3, then re-run this script.
""",
        encoding="utf-8",
    )

    print(f"Wrote {BUNDLE} ({len(manifest)} files copied)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
