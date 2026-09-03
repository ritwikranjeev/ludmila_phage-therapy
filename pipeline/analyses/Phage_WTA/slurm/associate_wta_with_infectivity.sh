#!/bin/bash
#SBATCH --job-name=WTA_ASSOC
#SBATCH --account=katrine_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

# pandas/sklearn via RCIC Anaconda; phage_wta_env prefix may be absent on compute nodes.
export PHAGE_WTA_SKIP_CONDA=1

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

# Toy / smoke: few hosts — use e.g.  export PHAGE_WTA_ASSOC_MIN_HOSTS=2  before sbatch, or pass flags by calling python3 directly.
# DRY_RUN=1: echo only (export DRY_RUN=1 before sbatch).
if [[ "${DRY_RUN:-0}" == "1" ]]; then
  echo "[DRY_RUN] Would run:"
  echo "  python3 \"${PHAGE_WTA_PROJECT_ROOT}/bin/associate_wta_with_infectivity.py\""
  echo "[DRY_RUN] Reads: data/wta_features/wta_feature_matrix.tsv, data/phenotypes/ludmila_infectivity.tsv"
  echo "[DRY_RUN] Writes: analysis/wta_infectivity_logreg_coeffs.tsv"
  echo "[DRY_RUN] Min merged hosts: \${PHAGE_WTA_ASSOC_MIN_HOSTS:-5} (set to 2 for toy tables)"
  exit 0
fi

python3 "${PHAGE_WTA_PROJECT_ROOT}/bin/associate_wta_with_infectivity.py"

echo "Done."
