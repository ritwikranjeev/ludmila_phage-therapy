#!/usr/bin/env bash
# Supplementary Tables 1-3 — the computed tables.
#
#   1  rate of adsorption, Welch by timepoint      ../fig04 data + helpers
#   2  one-step growth, Welch by timepoint         ../fig04 data + helpers
#   3  Kirby-Bauer, both zero policies             ../fig06 data + helpers
#
#   ./run.sh          all three
#   ./run.sh 3        one table
#
# These tables hold no data and no statistics of their own. Each calls the same
# function that produces the corresponding figure, so a table and its panel
# cannot report different numbers.
#
# Supplementary Tables 4 and 5 are authored by hand and are not generated here.
set -uo pipefail
cd "$(dirname "$0")"
mkdir -p logs output
pass=0; fail=0; want="${1:-}"
FIG04="$PWD/../fig04"
FIG06="$PWD/../fig06"

run_table() {
  local n="$1" script="$2" data="$3" shared="$4"
  [ -n "$want" ] && [ "$want" != "$n" ] && return
  if [ ! -e "$data" ] || [ ! -e "$shared" ]; then
    printf '  FAIL  %-1s %-26s missing input: %s\n' "$n" "$script" "$data"; fail=$((fail+1)); return
  fi
  if FIG_DATA="$data" FIG_OUT="$PWD/output" FIG_SHARED="$shared" \
     Rscript "scripts/$script" > "logs/table$n.log" 2>&1; then
    printf '  OK    %-1s %-26s -> output/\n' "$n" "$script"; pass=$((pass+1))
  else
    printf '  FAIL  %-1s %-26s -> logs/table%s.log\n' "$n" "$script" "$n"; fail=$((fail+1))
  fi
}

echo "Supplementary Tables"
run_table 1 supp_table01_roa.R    "$FIG04/data/4a" "$FIG04/scripts/_shared.R"
run_table 2 supp_table02_osgc.R   "$FIG04/data/4b" "$FIG04/scripts/_shared.R"
run_table 3 supp_table03_kirby.R  "$FIG06/data"    "$FIG06/scripts/_shared.R"
echo
echo "  passed $pass   failed $fail"
