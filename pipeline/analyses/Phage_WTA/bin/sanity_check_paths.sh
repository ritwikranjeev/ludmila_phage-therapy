#!/bin/bash
# Create key pipeline directories if missing (no Slurm jobs).
# Run from project root:  bash bin/sanity_check_paths.sh

set -euo pipefail

export PHAGE_WTA_SKIP_CONDA=1
# shellcheck source=/dev/null
source "$(dirname "$0")/../config/env_setup.sh"

ROOT="${PHAGE_WTA_PROJECT_ROOT}"
created=()

ensure_dir() {
  local d="$1"
  if [[ ! -d "$d" ]]; then
    mkdir -p "$d"
    created+=("$d")
  fi
}

ensure_dir "${ROOT}/logs"
ensure_dir "${ROOT}/data/phage"
ensure_dir "${ROOT}/data/hosts"
ensure_dir "${ROOT}/data/phage_prokka"
ensure_dir "${ROOT}/data/host_prokka"
ensure_dir "${ROOT}/data/phage_rbp"
ensure_dir "${ROOT}/data/wta_hmms"
ensure_dir "${ROOT}/data/wta_features/by_host"
ensure_dir "${ROOT}/data/alphafold/Ludmila_RBPs"
ensure_dir "${ROOT}/data/phenotypes"
ensure_dir "${ROOT}/analysis"

if [[ "${#created[@]}" -eq 0 ]]; then
  echo "All key directories already exist under ${ROOT}"
else
  echo "Created ${#created[@]} director(y/ies):"
  printf '  %s\n' "${created[@]}"
fi
