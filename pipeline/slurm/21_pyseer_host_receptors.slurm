#!/usr/bin/env bash
# Submit: sbatch slurm/21_pyseer_host_receptors.slurm
# Prereqs: 19 (host receptor Rtab), 13 (phenotypes), 07 (IQ-TREE tree), dendropy in pyseer env

#SBATCH -J pyseer_host_rec
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

bash "${PROJECT_ROOT}/scripts/21_pyseer_host_receptors.sh"
