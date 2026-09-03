#!/bin/bash
#SBATCH --job-name=WTA_HMMDB
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

# DRY_RUN=1: echo only (export DRY_RUN=1 before sbatch).
if [[ "${DRY_RUN:-0}" == "1" ]]; then
  echo "[DRY_RUN] Would run:"
  echo "  python3 \"${PHAGE_WTA_PROJECT_ROOT}/bin/build_wta_hmmdb.py\""
  echo "[DRY_RUN] Input: ${PHAGE_WTA_PROJECT_ROOT}/config/wta_reference_genes.fasta"
  echo "[DRY_RUN] Output dir: ${PHAGE_WTA_PROJECT_ROOT}/data/wta_hmms/ (per-gene *.hmm + wta_combined.hmm + hmmpress)"
  exit 0
fi

python3 "${PHAGE_WTA_PROJECT_ROOT}/bin/build_wta_hmmdb.py"

echo "Done."
