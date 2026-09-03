#!/usr/bin/env bash
# Supplementary Figure 4 — adsorption and burst size, pre- vs post-treatment isolate.
#
#   4a  adsorption rate   scripts/4a/adsorption_rate.R
#   4b  burst size        scripts/4b/burst_size.R
#
# Both panels share scripts/_shared.R.
#
#   ./run.sh          both panels
#   ./run.sh 4a       one panel
#
# Each panel writes, besides the figure:
#   *_welch.tsv          the headline comparison
#   *_window_sweep.tsv   the same comparison under every other defensible window
#
# Two constants were assumed rather than measured — host density for 4a and
# infective centres for 4b. Both cancel out of the pre/post ratio and both scale
# the absolute values, so they are printed in the log and recorded in the TSVs.
# To see the effect of a corrected value:
#
#   S4_HOST_DENSITY=5e7 ./run.sh 4a
#   S4_IC0=2e5 ./run.sh 4b
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

echo "Supplementary Figure 4"
run_panel 4a adsorption_rate.R
run_panel 4b burst_size.R
echo
echo "  passed $pass   failed $fail"
