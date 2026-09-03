#!/usr/bin/env bash
# Submit: sbatch slurm/05a_diagnose_contamination.slurm
# One task per sample in config/diagnose_contamination_samples.txt (default 3 samples → 0–2).

#SBATCH -J diag-contam
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 08:00:00
#SBATCH --array=0-2%2
#SBATCH -o outputs/logs/%x_%A_%a.out
#SBATCH -e outputs/logs/%x_%A_%a.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR}"
# shellcheck source=/dev/null
source "${PROJECT_ROOT}/config/hpc_env.sh"
hpc_setup_step_05a_diagnose_contamination

export THREADS="${SLURM_CPUS_PER_TASK:-8}"
bash "${PROJECT_ROOT}/scripts/05a_diagnose_contamination.sh"
