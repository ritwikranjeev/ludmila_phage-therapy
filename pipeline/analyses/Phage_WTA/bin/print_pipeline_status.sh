#!/bin/bash
# Print FOUND/MISSING for key Phage_WTA pipeline outputs (does not submit jobs).
# Run from project root:  bash bin/print_pipeline_status.sh

set -euo pipefail

# Status check only needs project root, not conda tools.
export PHAGE_WTA_SKIP_CONDA=1
# shellcheck source=/dev/null
source "$(dirname "$0")/../config/env_setup.sh"

ROOT="${PHAGE_WTA_PROJECT_ROOT}"

status_file() {
  local label="$1"
  local path="$2"
  if [[ -f "$path" ]]; then
    printf '%-55s FOUND   (%s)\n' "$label" "$path"
  else
    printf '%-55s MISSING (%s)\n' "$label" "$path"
  fi
}

status_dir_nonempty() {
  local label="$1"
  local path="$2"
  if [[ -d "$path" ]] && [[ -n "$(find "$path" -mindepth 1 -maxdepth 1 2>/dev/null | head -1)" ]]; then
    printf '%-55s FOUND   (%s; has entries)\n' "$label" "$path"
  else
    printf '%-55s MISSING (%s)\n' "$label" "$path"
  fi
}

echo "PHAGE_WTA_PROJECT_ROOT=${ROOT}"
echo "----------------------------------------"

status_file "config/phage_genomes.tsv" "${ROOT}/config/phage_genomes.tsv"
status_file "config/host_genomes.tsv" "${ROOT}/config/host_genomes.tsv"

status_dir_nonempty "data/host_prokka/ (subdirectories or files)" "${ROOT}/data/host_prokka"
status_dir_nonempty "data/phage_prokka/ (subdirectories or files)" "${ROOT}/data/phage_prokka"

status_file "data/wta_features/wta_feature_matrix.tsv" "${ROOT}/data/wta_features/wta_feature_matrix.tsv"
status_file "data/phenotypes/ludmila_infectivity.tsv" "${ROOT}/data/phenotypes/ludmila_infectivity.tsv"
status_file "analysis/wta_infectivity_logreg_coeffs.tsv" "${ROOT}/analysis/wta_infectivity_logreg_coeffs.tsv"

status_dir_nonempty "data/alphafold/Ludmila_RBPs" "${ROOT}/data/alphafold/Ludmila_RBPs"

echo "----------------------------------------"
echo "Done (no jobs submitted)."
