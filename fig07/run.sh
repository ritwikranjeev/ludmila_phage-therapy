#!/usr/bin/env bash
# Figure 7 — host range, and its supplementary detection-limit panel.
#
#   ./run.sh            draw the figure (R)
#   ./run.sh analysis   re-run the association analysis first (Python), then draw
#
# The figure is R, like the rest of the paper. The analysis behind it is Python
# because it needs pyhmmer for gene calling; it writes three small tables into
# work/ that the R script reads:
#
#   snp_distance.tsv     pairwise core-genome SNP distances
#   strain_metadata.tsv  group, susceptibility and source collection per strain
#   power_curve.tsv      simulated power by gene frequency and odds ratio
#
# Those tables are committed, so the figure redraws without re-running any of
# the heavy analysis. The analysis pipeline, in order, is:
#
#   scripts/recall_genes.py       gene calling      (~2 h, 8 threads)
#   scripts/association.py        mixed-model scan  (~1 min)
#   scripts/allele_analysis.py    WTA/LTA alleles   (~2 min)
set -uo pipefail
cd "$(dirname "$0")"
mkdir -p logs output work
HR="${HOST_RANGE:-$PWD/../../OCPT-Ludmila_Phage-Therapy/Host-Range}"

if [ "${1:-}" = "analysis" ]; then
  echo "Figure 7 — analysis"
  for step in "association.py" "allele_analysis.py" "gene_grid.py" "clonal_confounding.py" "supp_outputs.py"; do
    if python3 "scripts/$step" --host-range "$HR" --work work > "logs/${step%.py}.log" 2>&1; then
      printf '  OK    %-22s -> work/\n' "$step"
    else
      printf '  FAIL  %-22s -> logs/%s.log\n' "$step" "${step%.py}"
    fi
  done
fi

echo "Figure 7 — figure"
for f in work/gene_grid.tsv work/gene_grid_stats.tsv data/coregenome_tree.newick; do
  [ -f "$f" ] || { echo "  missing $f — run './run.sh analysis' first" >&2; exit 1; }
done

fail=0
# Figure 7 is a single figure; it has no supplementary panels. scripts/
# supp_figure_cluster_splitting.R is deliberately NOT drawn here. It is kept in
# the repository as documentation of why the gene re-calling step exists, and
# can be run by hand if a reviewer asks for it:
#
#   FIG_DATA=$PWD/data FIG_WORK=$PWD/work FIG_OUT=$PWD/output \
#     Rscript scripts/supp_figure_cluster_splitting.R
for r in figure7.R; do
  if FIG_DATA="$PWD/data" FIG_WORK="$PWD/work" FIG_OUT="$PWD/output" \
     Rscript "scripts/$r" > "logs/${r%.R}.log" 2>&1; then
    printf '  OK    %-34s -> output/\n' "$r"
  else
    printf '  FAIL  %-34s -> logs/%s.log\n' "$r" "${r%.R}"
    fail=1
  fi
done
exit $fail
