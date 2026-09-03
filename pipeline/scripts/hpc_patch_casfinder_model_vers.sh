#!/usr/bin/env bash
# One-time HPC fix: CasFinder models from `defense-finder update` may use vers="2.1"
# while MacSyFinder 2.1.4 (bundled with mdmparis-defense-finder) only accepts vers="2.0".
# Run on the login node after `defense-finder update` if you see:
#   MacsypyError: ... model definition ... has not the right version. version supported is '2.0'
#
# Usage: bash scripts/hpc_patch_casfinder_model_vers.sh
# Models dir defaults to ~/.macsyfinder/models/CasFinder/definitions

set -euo pipefail
ROOT="${MACSY_MODELS_CASFINDER:-${HOME}/.macsyfinder/models/CasFinder/definitions}"
[[ -d "${ROOT}" ]] || { echo "Missing ${ROOT}" >&2; exit 1; }
shopt -s nullglob
n=0
for f in "${ROOT}"/*.xml; do
  if grep -q 'vers="2.1"' "${f}" 2>/dev/null; then
    sed -i.bak 's/ vers="2.1"/ vers="2.0"/g' "${f}"
    n=$((n + 1))
  fi
done
echo "Patched vers 2.1->2.0 in ${n} XML files under ${ROOT}"
