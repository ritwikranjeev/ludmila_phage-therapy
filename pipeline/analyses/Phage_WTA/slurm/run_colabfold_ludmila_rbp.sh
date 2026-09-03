#!/bin/bash
#SBATCH --job-name=LUDMILA_RBP_AF
# HPC3 / RCIC: using free-gpu — no dedicated GPU Slurm account required (see RCIC Slurm docs).
# Jobs may be PREEMPTED; checkpoint/restart ColabFold if possible. For non-preempt runs, switch to
# partition=gpu and set --account to a lab GPU account (name often ends in _gpu) once your PI/RCIC adds one.
#SBATCH --account=katrine_lab
#SBATCH --partition=free-gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
# RCIC examples: --gres=gpu:1 or --gres=gpu:V100:1 (see `sinfo` / RCIC Slurm GPU examples).
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
# Max walltime for free-gpu may differ from gpu (see: scontrol show partition free-gpu).
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

# --- GPU / CUDA (HPC3): conda ColabFold often bundles CUDA; only add modules if needed. ---
# Check: module spider cuda    module spider cudnn
# Examples (uncomment after verifying versions on HPC3):
# module load cuda/12.x
# module load cudnn/9.x

INPUT_FASTA="${PHAGE_WTA_PROJECT_ROOT}/data/phage_rbp/Ludmila_RBP_candidates.fasta"
OUTPUT_DIR="${PHAGE_WTA_PROJECT_ROOT}/data/alphafold/Ludmila_RBPs"

mkdir -p "${OUTPUT_DIR}"

hostname
date
echo "INPUT_FASTA=${INPUT_FASTA}"
echo "OUTPUT_DIR=${OUTPUT_DIR}"

# DRY_RUN=1: echo only (export DRY_RUN=1 before sbatch). Does not require INPUT_FASTA to exist.
if [[ "${DRY_RUN:-0}" == "1" ]]; then
  echo "[DRY_RUN] Would run (add flags per ColabFold docs / TODO in this script):"
  echo "  colabfold_batch \"${INPUT_FASTA}\" \"${OUTPUT_DIR}\""
  if [[ ! -f "${INPUT_FASTA}" ]]; then
    echo "[DRY_RUN] NOTE: INPUT_FASTA missing — create with extract_ludmila_rbp_candidates.py before real run."
  fi
  exit 0
fi

if [[ ! -f "${INPUT_FASTA}" ]]; then
  echo "ERROR: Missing candidate FASTA: ${INPUT_FASTA}" >&2
  echo "Run: python3 bin/extract_ludmila_rbp_candidates.py ..." >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# ColabFold binary: install in an env visible from GPU nodes, or set COLABFOLD_BATCH
# to the full path (e.g. under /dfs6b/... if login-only conda prefix breaks on workers).
# Optional: uncomment CUDA/cudnn modules above if the bundled runtime is insufficient.
# Add flags as needed: --model-type, --num-recycle, etc. (see ColabFold docs).
# ---------------------------------------------------------------------------
CFB_CMD=(colabfold_batch)
if [[ -n "${COLABFOLD_BATCH:-}" ]]; then
  CFB_CMD=("${COLABFOLD_BATCH}")
elif ! command -v colabfold_batch &>/dev/null; then
  echo "ERROR: colabfold_batch not on PATH and COLABFOLD_BATCH is unset." >&2
  echo "  Install ColabFold (conda/pip per upstream) on a filesystem mounted on GPU nodes, then:" >&2
  echo "    export COLABFOLD_BATCH=/path/to/colabfold_batch" >&2
  echo "  Or ensure phage_wta_env (or another env) is activated and provides colabfold_batch." >&2
  exit 1
fi

"${CFB_CMD[@]}" "${INPUT_FASTA}" "${OUTPUT_DIR}"

echo "Done."
