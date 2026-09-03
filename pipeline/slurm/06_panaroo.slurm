#!/usr/bin/env bash
# Submit: sbatch slurm/06_panaroo.slurm
#SBATCH -J panaroo
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -o outputs/logs/%x_%A.out
#SBATCH -e outputs/logs/%x_%A.err
set -euo pipefail
PROJECT_ROOT="${SLURM_SUBMIT_DIR}"
# shellcheck source=/dev/null
source "${PROJECT_ROOT}/config/hpc_env.sh"
hpc_setup_step_06_panaroo

bash "${PROJECT_ROOT}/scripts/06_panaroo.sh"
