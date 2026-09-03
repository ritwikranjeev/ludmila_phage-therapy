#!/usr/bin/env bash
# Submit: sbatch slurm/04d_quast_cleaned.slurm

#SBATCH -J quast-cln
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -o outputs/logs/%x_%A.out
#SBATCH -e outputs/logs/%x_%A.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR}"
# shellcheck source=/dev/null
source "${PROJECT_ROOT}/config/hpc_env.sh"
hpc_setup_step_04_quast

bash "${PROJECT_ROOT}/scripts/04d_quast_cleaned.sh"
