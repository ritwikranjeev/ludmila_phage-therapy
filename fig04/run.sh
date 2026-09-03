#!/usr/bin/env bash
# Figure 4 — reduced efficiency of OCPT-Ludmila against the post-treatment isolate.
#
#   4a  rate of attachment      scripts/4a/roa_attached.R
#   4b  one-step growth         scripts/4b/one_step_growth.R
#   4c  efficiency of plating   scripts/4c/plot_eop_pre_post.R
#
# Panels 4a and 4b share plotting helpers via scripts/_shared.R. This file is the
# only place that knows the directory layout.
#
#   ./run.sh          all panels
#   ./run.sh 4b       one panel
set -uo pipefail
cd "$(dirname "$0")"
mkdir -p logs
pass=0; fail=0; want="${1:-}"

run_panel() {
  local panel="$1" script="$2"
  [ -n "$want" ] && [ "$want" != "$panel" ] && return
  if FIG_DATA="$PWD/data/$panel" FIG_OUT="$PWD/output/$panel" \
     FIG_SHARED="$PWD/scripts/_shared.R" \
     Rscript "scripts/$panel/$script" > "logs/$panel.log" 2>&1; then
    printf '  OK    %-3s %-22s -> output/%s/\n' "$panel" "$script" "$panel"; pass=$((pass+1))
  else
    printf '  FAIL  %-3s %-22s -> logs/%s.log\n' "$panel" "$script" "$panel"; fail=$((fail+1))
  fi
}

echo "Figure 4"
run_panel 4a roa_attached.R
run_panel 4b one_step_growth.R
run_panel 4c plot_eop_pre_post.R
echo
echo "  passed $pass   failed $fail"
