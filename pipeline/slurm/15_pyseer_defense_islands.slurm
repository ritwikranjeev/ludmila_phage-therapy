#!/usr/bin/env bash
# Submit: sbatch slurm/15_pyseer_defense_islands.slurm
#SBATCH -J pyseer_defense
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 08:00:00
#SBATCH -o outputs/logs/%x_%A.out
#SBATCH -e outputs/logs/%x_%A.err
set -euo pipefail
PROJECT_ROOT="${SLURM_SUBMIT_DIR}"
# shellcheck source=/dev/null
source "${PROJECT_ROOT}/config/hpc_env.sh"
hpc_setup_step_pyseer

bash "${PROJECT_ROOT}/scripts/15_pyseer_defense_islands.sh"
