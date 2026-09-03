#!/usr/bin/env bash
# Supplementary Figure 5 — mutation maps for the remaining host backgrounds,
# and the two hypermutator replicates shown on their own.
#
#   5a  USA300 mutation map                 FIG_STRAINS=USA300      FIG_MAPS=main,shaped
#   5b  SA412 mutation map                  FIG_STRAINS=SA412       FIG_MAPS=main,shaped
#   5c  PedCF mutation map                  FIG_STRAINS=PedCF       FIG_MAPS=main,shaped
#   5d  OCPT04SA01 hypermutator (rep. 3)    FIG_STRAINS=OCPT04SA01  FIG_MAPS=hypermutator
#   5e  PedCF hypermutator (rep. 2)         FIG_STRAINS=PedCF       FIG_MAPS=hypermutator
#
# Generator and inputs live in ../_shared/strain_maps: the same breseq output
#
# The breseq / Clair3 workflows that produced those mutation calls are in
# ../_shared/variant_calling/ (in_vitro_mutants for Illumina, clinical_ont for ONT).
# feeds Figure 5a/5c and every panel here, so there is one copy of it. This
# figure directory holds no data of its own by design.
#
# The two hypermutator replicates are excluded from the multi-replicate maps
# (their variant counts would swamp the other replicates) and drawn here
# instead. That exclusion is declared in _strain_maps.py as HYPERMUTATOR_SAMPLES.
#
# UniProt lookups are OFF by default and served from the committed cache.
# Set FIG_ALLOW_NETWORK=1 only to refresh it.
#
#   ./run.sh          all panels
#   ./run.sh 5d       one panel
set -uo pipefail
cd "$(dirname "$0")"
mkdir -p logs
pass=0; fail=0; want="${1:-}"
PY=${PYTHON:-python3}
SHARED="$PWD/../_shared/strain_maps"

panel() {
  local panel="$1" strain="$2" maps="$3"
  [ -n "$want" ] && [ "$want" != "$panel" ] && return
  if FIG_DATA="$SHARED/data" FIG_OUT="$PWD/output/$panel" \
     FIG_STRAINS="$strain" FIG_MAPS="$maps" \
     "$PY" "$SHARED/_strain_maps.py" > "logs/$panel.log" 2>&1; then
    printf '  OK    %-3s %-32s -> output/%s/\n' "$panel" "$strain ($maps)" "$panel"; pass=$((pass+1))
  else
    printf '  FAIL  %-3s %-32s -> logs/%s.log\n' "$panel" "$strain ($maps)" "$panel"; fail=$((fail+1))
  fi
}

echo "Supplementary Figure 5"
panel 5a USA300     main,shaped
panel 5b SA412      main,shaped
panel 5c PedCF      main,shaped
panel 5d OCPT04SA01 hypermutator
panel 5e PedCF      hypermutator
echo
echo "  passed $pass   failed $fail"
