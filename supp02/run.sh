#!/usr/bin/env bash
# Supplementary Figure 2 — genomics and stability.
#
#   2a  proteomic tree (ViPTree)     scripts/2a/proteomic_tree.R
#   2b  pH stability                 scripts/2b/ph_stability.R
#   2c  thermal stability            scripts/2c/thermal_stability.R
#   2d  VIRIDIC similarity heatmap   scripts/2d/viridic_heatmap.R
#
# Panels 2b and 2c share scripts/_shared.R. Panels 2a and 2d are standalone and
# carry heavier dependencies: 2a needs ape and ggtree (Bioconductor), 2d needs
# only ggplot2.
#
#   ./run.sh          all panels
#   ./run.sh 2c       one panel
#
# 2a renders both a rectangular and a circular layout, since which one goes in
# the supplement is a layout choice rather than an analysis:
#
#   S2_LAYOUTS=rectangular ./run.sh 2a
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
    printf '  OK    %-3s %-24s -> output/%s/\n' "$panel" "$script" "$panel"; pass=$((pass+1))
  else
    printf '  FAIL  %-3s %-24s -> logs/%s.log\n' "$panel" "$script" "$panel"; fail=$((fail+1))
  fi
}

echo "Supplementary Figure 2"
run_panel 2a proteomic_tree.R
run_panel 2b ph_stability.R
run_panel 2c thermal_stability.R
run_panel 2d viridic_heatmap.R
echo
echo "  passed $pass   failed $fail"
