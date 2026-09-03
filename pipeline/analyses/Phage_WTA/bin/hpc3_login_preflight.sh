#!/usr/bin/env bash
# Run preflight on an HPC3 login node when default `python3` is too old (e.g. 3.6 from ~/.local mamba).
# Loads an RCIC Anaconda module for a modern interpreter, then runs preflight with PHAGE_WTA_SKIP_CONDA=1.
#
#   cd /pub/ritwik1/manuscripts/OCPT-Ludmila_PhageTherapy/WTA-RBP_Analysis
#   bash bin/hpc3_login_preflight.sh --allow-wta-placeholder

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "${ROOT}"

if command -v module &>/dev/null; then
  _loaded=
  for _mod in anaconda/2024.06 anaconda/2025.12 anaconda/2022.05 anaconda/2021.11; do
    if module load "${_mod}" 2>/dev/null; then
      _loaded=1
      break
    fi
  done
  unset _mod
  if [[ -z "${_loaded}" ]]; then
    echo "WARNING: Could not module load anaconda (tried 2024.06, 2025.12, 2022.05, 2021.11). Set PHAGE_WTA_CONDA_MODULE." >&2
  fi
  unset _loaded
fi

export PHAGE_WTA_SKIP_CONDA=1
# shellcheck source=/dev/null
source "${ROOT}/config/env_setup.sh"

exec python3 "${ROOT}/bin/preflight_check.py" "$@"
