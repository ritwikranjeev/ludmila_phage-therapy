#!/usr/bin/env bash
# Figure 2 — morphological, genomic and kinetic characterisation of OCPT-Ludmila.
#
# This file is the only place that knows the directory layout. Each script is
# handed its input and output directory and knows nothing else.
#
#   ./run.sh            run every panel
#   ./run.sh 2d         run one panel
set -uo pipefail
cd "$(dirname "$0")"
mkdir -p logs

# panel : script
declare -a PANELS=(
  "2d:osgc.R"
  "2e:moi_curves.R"
)
# 2a and 2b are micrographs with no generating code.
#
# 2c is not regenerated here either. The deposited genome map,
# output/2c/ludmila_phold_circos.png and .svg, was written by Phold's own
# plotting routine as step 07 of scripts/2c/genome_pipeline/. Reproducing it
# needs the Phold databases and a structure-prediction pass, so it cannot run
# from this repository; the pipeline scripts are deposited as the provenance
# record instead.

pass=0; fail=0; skip=0
run_panel() {
  local panel="$1" script="$2"
  local sp="scripts/$panel/$script"
  if [ ! -f "$sp" ]; then printf '  SKIP  %-4s %-24s (no script)\n' "$panel" "$script"; skip=$((skip+1)); return; fi
  if [ ! -d "data/$panel" ] || [ -z "$(ls -A "data/$panel" 2>/dev/null)" ]; then
    printf '  SKIP  %-4s %-24s (data/%s is empty)\n' "$panel" "$script" "$panel"; skip=$((skip+1)); return; fi
  if FIG_DATA="$PWD/data/$panel" FIG_OUT="$PWD/output/$panel" \
     Rscript "$sp" > "logs/$panel.log" 2>&1; then
    printf '  OK    %-4s %-24s -> output/%s/\n' "$panel" "$script" "$panel"; pass=$((pass+1))
  else
    printf '  FAIL  %-4s %-24s -> logs/%s.log\n' "$panel" "$script" "$panel"; fail=$((fail+1))
  fi
}

# --- optional: confirm 2e still matches the original hardcoded script ---------
if [ "${1:-}" = "--verify-2e" ]; then
  FIG_DATA="$PWD/data/2e" Rscript scripts/2e/verify_against_original.R
  exit $?
fi

want="${1:-}"
echo "Figure 2"
for entry in "${PANELS[@]}"; do
  panel="${entry%%:*}"; script="${entry#*:}"
  [ -n "$want" ] && [ "$want" != "$panel" ] && continue
  run_panel "$panel" "$script"
done
echo
echo "  passed $pass   failed $fail   skipped $skip"
