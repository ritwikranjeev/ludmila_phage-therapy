#!/usr/bin/env bash
# Submit: sbatch slurm/05_prokka.slurm
# Array 0-82, max 20 concurrent (%20).
#SBATCH -J prokka
#SBATCH --array=0-82%20
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH -o outputs/logs/%x_%A_%a.out
#SBATCH -e outputs/logs/%x_%A_%a.err
set -euo pipefail
PROJECT_ROOT="${SLURM_SUBMIT_DIR}"
# shellcheck source=/dev/null
source "${PROJECT_ROOT}/config/hpc_env.sh"
hpc_setup_step_05_prokka

bash "${PROJECT_ROOT}/scripts/05_prokka.sh"
