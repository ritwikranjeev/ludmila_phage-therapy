#!/bin/bash
# Copy this file for new GPU jobs. DRY_RUN=1 (export before sbatch) prints example colabfold_batch and exits 0.
# Default partition free-gpu matches slurm/run_colabfold_ludmila_rbp.sh (no GPU lab account). Use gpu + *_gpu account when available.
#SBATCH --job-name=TEST_GPU
#SBATCH --account=katrine_lab
#SBATCH --partition=free-gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=08:00:00
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

# --- GPU / CUDA modules (HPC3): fill in exact names from RCIC documentation / `module spider cuda` ---
# module load cuda/...
# module load cudnn/...

hostname
date

INPUT_FASTA="${PHAGE_WTA_PROJECT_ROOT}/data/phage_rbp/Ludmila_RBP_candidates.fasta"
OUTPUT_DIR="${PHAGE_WTA_PROJECT_ROOT}/data/alphafold/Ludmila_RBPs"

if [[ "${DRY_RUN:-0}" == "1" ]]; then
  echo "[DRY_RUN] Example GPU workload (ColabFold):"
  echo "  colabfold_batch \"${INPUT_FASTA}\" \"${OUTPUT_DIR}\""
  echo "[DRY_RUN] Uncomment module load cuda/cudnn above if your install needs cluster modules."
  exit 0
fi

echo "Replace with colabfold_batch input.fasta output_dir"
