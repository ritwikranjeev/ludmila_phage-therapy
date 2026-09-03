#!/usr/bin/env bash
# Supplementary Figure 3 — serum neutralisation. One script, one output, no panels.
#
#   ./run.sh
#
# EOP is computed from plaque counts relative to pre-treatment serum (week 0).
# The stats table reports both the Welch test the published analysis used and a
# count-based Poisson test, plus the log10 reduction; they disagree, and the
# log reduction is the scale neutralisation is judged on.
set -uo pipefail
cd "$(dirname "$0")"
mkdir -p logs

echo "Supplementary Figure 3"
if FIG_DATA="$PWD/data" FIG_OUT="$PWD/output" \
   Rscript scripts/serum_neutralisation.R > logs/supp03.log 2>&1; then
  printf '  OK    serum_neutralisation.R -> output/\n'
else
  printf '  FAIL  serum_neutralisation.R -> logs/supp03.log\n'
  exit 1
fi
