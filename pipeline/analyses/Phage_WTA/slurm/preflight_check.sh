#!/bin/bash
#SBATCH --job-name=WTA_PREFLT
#SBATCH --account=katrine_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

if [[ "${DRY_RUN:-0}" == "1" ]]; then
  export PHAGE_WTA_SKIP_CONDA=1
fi

# Slurm copies the batch script to the spool directory; paths from $0 break. Submit from project root.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  cd "${SLURM_SUBMIT_DIR}"
else
  cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
# shellcheck source=/dev/null
source "${PWD}/config/env_setup.sh"

cd "${PHAGE_WTA_PROJECT_ROOT}"
mkdir -p logs

hostname
date

# On cluster, fail if WTA still has placeholders unless you pass this export:
#   sbatch --export=ALL,PREFLIGHT_ALLOW_WTA_PLACEHOLDER=1 slurm/preflight_check.sh
EXTRA=()
if [[ "${PREFLIGHT_ALLOW_WTA_PLACEHOLDER:-0}" == "1" ]]; then
  EXTRA+=(--allow-wta-placeholder)
fi

# DRY_RUN=1: echo only (export DRY_RUN=1 before sbatch so it reaches the compute node).
if [[ "${DRY_RUN:-0}" == "1" ]]; then
  echo "[DRY_RUN] Would run:"
  if [[ "${#EXTRA[@]}" -eq 0 ]]; then
    echo "  python3 \"${PHAGE_WTA_PROJECT_ROOT}/bin/preflight_check.py\""
  else
    echo "  python3 \"${PHAGE_WTA_PROJECT_ROOT}/bin/preflight_check.py\" ${EXTRA[*]}"
  fi
  echo "[DRY_RUN] Set PREFLIGHT_ALLOW_WTA_PLACEHOLDER=1 to add --allow-wta-placeholder on real run."
  exit 0
fi

python3 "${PHAGE_WTA_PROJECT_ROOT}/bin/preflight_check.py" "${EXTRA[@]}"
