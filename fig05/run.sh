#!/usr/bin/env bash
# Figure 5 — convergent resistance evolution, in vitro and in vivo.
#
#   5a  OCPT04SA01 mutation map   ../_shared/strain_maps  (FIG_STRAINS=OCPT04SA01)
#   5b  clinical pre/post map     scripts/5b/clinical_pre_post_map.py
#   5c  OCPT02SA01 mutation map   ../_shared/strain_maps  (FIG_STRAINS=OCPT02SA01)
#
# The generator and its inputs live in ../_shared/strain_maps because the same
#
# The breseq / Clair3 workflows that produced those mutation calls are in
# ../_shared/variant_calling/ (in_vitro_mutants for Illumina, clinical_ont for ONT).
# breseq output feeds Figure 5a/5c and Supplementary Figure 5a-e. One copy, so
# the two figures cannot drift apart. FIG_MAPS picks which map a panel wants.
#
# This figure is still Python. 5b reuses the drawing machinery in _strain_maps.py,
# so that file is both the 5a/5c generator and 5b's shared library; it is passed
# to 5b via FIG_SHARED. Porting to R is on the outstanding list.
#
# UniProt lookups are OFF by default and served from the committed cache.
# Set FIG_ALLOW_NETWORK=1 only to refresh it.
#
#   ./run.sh          all panels
#   ./run.sh 5b       one panel
set -uo pipefail
cd "$(dirname "$0")"
mkdir -p logs
pass=0; fail=0; want="${1:-}"
PY=${PYTHON:-python3}
SHARED="$PWD/../_shared/strain_maps"

strain_panel() {
  local panel="$1" strain="$2"
  [ -n "$want" ] && [ "$want" != "$panel" ] && return
  if FIG_DATA="$SHARED/data" FIG_OUT="$PWD/output/$panel" FIG_STRAINS="$strain" \
     FIG_MAPS="main,shaped" \
     "$PY" "$SHARED/_strain_maps.py" > "logs/$panel.log" 2>&1; then
    printf '  OK    %-3s %-28s -> output/%s/\n' "$panel" "_strain_maps.py ($strain)" "$panel"; pass=$((pass+1))
  else
    printf '  FAIL  %-3s %-28s -> logs/%s.log\n' "$panel" "_strain_maps.py ($strain)" "$panel"; fail=$((fail+1))
  fi
}

echo "Figure 5"
strain_panel 5a OCPT04SA01

if [ -z "$want" ] || [ "$want" = "5b" ]; then
  if FIG_DATA="$PWD/data/5b" FIG_OUT="$PWD/output/5b" \
     FIG_SHARED="$SHARED/_strain_maps.py" \
     "$PY" scripts/5b/clinical_pre_post_map.py > logs/5b.log 2>&1; then
    printf '  OK    5b  %-28s -> output/5b/\n' "clinical_pre_post_map.py"; pass=$((pass+1))
  else
    printf '  FAIL  5b  %-28s -> logs/5b.log\n' "clinical_pre_post_map.py"; fail=$((fail+1))
  fi
fi

strain_panel 5c OCPT02SA01
echo
echo "  passed $pass   failed $fail"
