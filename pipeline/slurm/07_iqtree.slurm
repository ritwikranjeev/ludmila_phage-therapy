#!/usr/bin/env bash
# Submit: sbatch slurm/07_iqtree.slurm
#SBATCH -J iqtree
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH -o outputs/logs/%x_%A.out
#SBATCH -e outputs/logs/%x_%A.err
set -euo pipefail
PROJECT_ROOT="${SLURM_SUBMIT_DIR}"
# shellcheck source=/dev/null
source "${PROJECT_ROOT}/config/hpc_env.sh"
hpc_setup_step_07_iqtree

bash "${PROJECT_ROOT}/scripts/07_iqtree.sh"
