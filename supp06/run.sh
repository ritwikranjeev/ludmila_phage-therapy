#!/usr/bin/env bash
# Supplementary Figure 6 — Kirby-Bauer, all strains in one 3 x 2 grid.
#
#   ./run.sh
#
# This figure holds no data and no statistics of its own. It reads Figure 6's
# inputs and Figure 6's helpers, so the test and the correction stay identical
# between the main figure and the supplement.
#
#   ../fig06/data/6a/kirby_bauer_in_vitro_mutants.tsv
#   ../fig06/data/6b/kirby_bauer_pre_post.tsv
#   ../fig06/scripts/_shared.R
#
# Benjamini-Hochberg is applied across all 36 comparisons in the grid, not
# within each panel: the reader sees them as one display item.
#
# Same switches as Figure 6, e.g.:
#   KB_TEST=wilcoxon ./run.sh            # rank-based instead of Welch
set -uo pipefail
cd "$(dirname "$0")"
mkdir -p logs
FIG06="$PWD/../fig06"

for f in "$FIG06/data/6a" "$FIG06/data/6b" "$FIG06/scripts/_shared.R"; do
  if [ ! -e "$f" ]; then
    echo "  Missing Figure 6 input: $f" >&2
    exit 1
  fi
done

echo "Supplementary Figure 6"
if FIG_DATA="$FIG06/data" FIG_OUT="$PWD/output" \
   FIG_SHARED="$FIG06/scripts/_shared.R" \
   Rscript scripts/kirby_combined.R > logs/supp06.log 2>&1; then
  printf '  OK    kirby_combined.R -> output/\n'
else
  printf '  FAIL  kirby_combined.R -> logs/supp06.log\n'
  exit 1
fi
