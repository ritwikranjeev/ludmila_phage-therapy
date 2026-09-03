#!/usr/bin/env bash
# Submit: sbatch slurm/01_fetch_reads.slurm
# Array 0-82, max 20 concurrent tasks (%20). Override: sbatch --array=0-82%N ...

#SBATCH -J fetch_reads
#SBATCH --array=0-82%20
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 04:00:00
#SBATCH -o outputs/logs/%x_%A_%a.out
#SBATCH -e outputs/logs/%x_%A_%a.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR}"
# shellcheck source=/dev/null
source "${PROJECT_ROOT}/config/hpc_env.sh"
hpc_setup_step_01_fetch

bash "${PROJECT_ROOT}/scripts/01_fetch_reads.sh"
