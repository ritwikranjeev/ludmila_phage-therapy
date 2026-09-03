#!/usr/bin/env python3
"""Fill perplexity_upload_batch/ with copies for batch upload + paste-ready prompt."""

import shutil
from pathlib import Path
from datetime import datetime, timezone

PROJECT_ROOT = Path(__file__).resolve().parents[1]
BUNDLE = PROJECT_ROOT / "outputs" / "analysis_planning_bundle"
OUT = PROJECT_ROOT / "perplexity_upload_batch"

PROMPT = """CONTEXT FOR FOLLOW-UP DISCUSSION (paste this message into Perplexity, then attach all files in the folder)

I am planning the next phase of analysis for a bacterial genome-wide association study (GWAS) project. Attached are tab-delimited tables and short READMEs from our pipeline.

Study setup:
- Host: Staphylococcus aureus isolates (see samples.tsv; n ≈ 83) with phage infection phenotypes (binary and 0–3 scale; see host_range.tsv and derived phenotypes.tsv).
- Association: susceptibility vs resistance to a podovirus phage (“Ludmila”). Phage genome annotated with phold; tail/receptor-associated CDS flagged in phage_receptor_associated.tsv.
- GWAS: pyseer with patristic distances from a core-genome IQ-TREE tree (population_structure_dist.tsv) to control structure. Separate runs: pangenome genes, DefenseFinder families, prophage binary features (PhiSpy), curated WTA / surface receptor panel (host_receptor_gene_panel.tsv).
- FDR: Benjamini–Hochberg q-values were computed per pyseer run (not one global merge across all layers). gwas_priority_overview.tsv is an exploratory merge of top rows only.

Please help me with:
1) How to prioritize experimental or computational follow-up on host hits (genes, defense, WTA loci) given multiple testing and unstable pyseer “notes” rows.
2) How to frame phage tail / RBP candidates vs host WTA diversity for a microbiology / infection audience (hypothesis text + figure ideas).
3) Caveats specific to sparse binary features (e.g. prophage counts) and defense GWAS.
4) A concise outline for main-text vs supplementary tables/figures.

If anything in the tables is unclear, say what extra metadata you would need.
"""


def main() -> int:
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)

    copies = []

    def cpy(rel_src, dest_name=None):
        if not rel_src.is_file():
            return
        name = dest_name or rel_src.name
        dst = OUT / name
        shutil.copy2(rel_src, dst)
        copies.append(name)

    # From bundle (preferred)
    if BUNDLE.is_dir():
        for sub in (
            "config_and_phenotypes/samples.tsv",
            "config_and_phenotypes/phenotypes.tsv",
            "config_and_phenotypes/host_range.tsv",
            "config_and_phenotypes/host_receptor_gene_panel.tsv",
            "gwas/gwas_priority_overview.tsv",
            "gwas/population_structure_dist.tsv",
            "integration/phage_host_receptor_crosswalk.tsv",
            "integration/host_receptor_panel_reference.tsv",
            "integration/host_receptor_gwas_top_binary.tsv",
            "phage/phage_genes.tsv",
            "phage/phage_receptor_associated.tsv",
            "prophages/prophage_summary.tsv",
            "prophages/prophage_feature_matrix.tsv",
            "defense/defense_islands.tsv",
            "defense/defense_family_presence.tsv",
            "trees/core.treefile",
        ):
            cpy(BUNDLE / sub)
        if (BUNDLE / "README.txt").is_file():
            cpy(BUNDLE / "README.txt", "00_analysis_planning_bundle_README.txt")
        if (BUNDLE / "MANIFEST.txt").is_file():
            cpy(BUNDLE / "MANIFEST.txt", "00_MANIFEST_file_sources.txt")
        hits = BUNDLE / "gwas" / "top_hits"
        if hits.is_dir():
            for p in sorted(hits.glob("*.tsv")):
                cpy(p, f"gwas_top_{p.name}")
        fdr = BUNDLE / "gwas" / "fdr_tables"
        if fdr.is_dir():
            subdir = OUT / "gwas_fdr_tables"
            subdir.mkdir()
            for p in sorted(fdr.glob("*.tsv")):
                shutil.copy2(p, subdir / p.name)
                copies.append(f"gwas_fdr_tables/{p.name}")
    else:
        # Minimal fallback from repo paths
        for rel in (
            "config/samples.tsv",
            "config/phenotypes.tsv",
            "data/phenotypes/host_range.tsv",
            "config/host_receptor_gene_panel.tsv",
        ):
            cpy(PROJECT_ROOT / rel)

    (OUT / "START_HERE_PASTE_INTO_PERPLEXITY.txt").write_text(PROMPT.strip() + "\n", encoding="utf-8")
    copies.insert(0, "START_HERE_PASTE_INTO_PERPLEXITY.txt")

    idx = [
        f"# perplexity_upload_batch — built {datetime.now(timezone.utc).isoformat()}",
        "# Upload: select all files in this folder into Perplexity, or zip from repo root:",
        "#   zip -r perplexity_upload_batch.zip perplexity_upload_batch/",
        "# Paste the text from START_HERE_PASTE_INTO_PERPLEXITY.txt as your first message.",
        "",
        "Files:",
    ]
    for n in sorted(set(copies)):
        idx.append(f"  - {n}")
    if (OUT / "gwas_fdr_tables").exists():
        for p in sorted((OUT / "gwas_fdr_tables").glob("*.tsv")):
            idx.append(f"  - gwas_fdr_tables/{p.name}")

    (OUT / "00_FILE_INDEX.txt").write_text("\n".join(idx) + "\n", encoding="utf-8")

    print(f"Wrote {OUT} ({len(list(OUT.rglob('*')))} paths)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
