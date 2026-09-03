#!/usr/bin/env bash
# Figure 6 — antibiotic susceptibility of phage-resistant mutants and of the
# pre/post clinical isolate pair.
#
#   6a  in vitro resistant mutants vs ancestors   scripts/6a/kirby_in_vitro_mutants.R
#   6b  pre- vs post-treatment clinical isolates  scripts/6b/kirby_pre_post.R
#
# Both panels share scripts/_shared.R. This file is the only place that knows
# the directory layout.
#
#   ./run.sh          both panels
#   ./run.sh 6a       one panel
#
# Zone diameters follow CLSI: the diameter includes the 6 mm disc, so a plate
# with no visible clearing reads 6 mm and no reading can fall below it. The
# scripts stop if one does, since a sub-6 mm value means the sheet is on a
# different scale. KB_DISC_MM overrides the floor if a different disc is used.
#
# An earlier version of this sheet recorded no zone as 0 mm, which made "no
# zone" and "not measured" indistinguishable; the published panels deleted those
# readings, which is what produced the two significant in vitro results. The
# re-measured sheet removes the ambiguity, so no zero policy is needed.
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
    printf '  OK    %-3s %-28s -> output/%s/\n' "$panel" "$script" "$panel"; pass=$((pass+1))
  else
    printf '  FAIL  %-3s %-28s -> logs/%s.log\n' "$panel" "$script" "$panel"; fail=$((fail+1))
  fi
}

echo "Figure 6"
run_panel 6a kirby_in_vitro_mutants.R
run_panel 6b kirby_pre_post.R
echo
echo "  passed $pass   failed $fail"
