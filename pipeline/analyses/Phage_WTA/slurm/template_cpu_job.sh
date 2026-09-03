#!/bin/bash
# Copy this file for new CPU jobs. DRY_RUN=1 (export before sbatch) prints example commands and exits 0.
#SBATCH --job-name=TEST_CPU
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

if [[ "${DRY_RUN:-0}" == "1" ]]; then
  echo "[DRY_RUN] Example CPU workload (replace with your command):"
  echo "  python3 \"\${PHAGE_WTA_PROJECT_ROOT}/bin/your_script.py\""
  echo "  # or: prokka --cpus \"\${SLURM_CPUS_PER_TASK:-8}\" --outdir ... --prefix ID genome.fa"
  exit 0
fi

echo "Replace this with a real command"
