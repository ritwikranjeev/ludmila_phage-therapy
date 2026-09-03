#!/usr/bin/env bash
# Figure 3 — treatment timeline and longitudinal laboratory values.
#
# All eight panels (a-h) come from one script as a single composite, so there
# are no per-panel directories. This file is the only place that knows the layout.
#
#   ./run.sh          regenerate the figure
set -uo pipefail
cd "$(dirname "$0")"
mkdir -p logs

if FIG_DATA="$PWD/data" FIG_OUT="$PWD/output" \
   Rscript scripts/treatment_timeline_figure.R > logs/fig03.log 2>&1; then
  echo "  OK    treatment_timeline_figure.R -> output/"
  ls output | sed 's/^/          /'
else
  echo "  FAIL  treatment_timeline_figure.R -> logs/fig03.log"
  tail -5 logs/fig03.log | sed 's/^/          /'
  exit 1
fi
