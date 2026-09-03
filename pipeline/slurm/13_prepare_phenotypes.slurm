#!/usr/bin/env bash
# Submit: sbatch slurm/13_prepare_phenotypes.slurm
#SBATCH -J prepare_pheno
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH -o outputs/logs/%x_%A.out
#SBATCH -e outputs/logs/%x_%A.err
set -euo pipefail
PROJECT_ROOT="${SLURM_SUBMIT_DIR}"
# shellcheck source=/dev/null
source "${PROJECT_ROOT}/config/hpc_env.sh"
hpc_setup_step_python_stdlib

python3 "${PROJECT_ROOT}/scripts/13_prepare_phenotypes.py"
