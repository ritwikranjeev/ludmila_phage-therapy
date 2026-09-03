#!/usr/bin/env bash
# Submit: sbatch slurm/04c_clean_assemblies.slurm
# One array task per sample in config/assembly_cleanup_samples.txt (currently 11 lines → 0–10).

#SBATCH -J asm-clean
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 04:00:00
#SBATCH --array=0-10%4
#SBATCH -o outputs/logs/%x_%A_%a.out
#SBATCH -e outputs/logs/%x_%A_%a.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR}"
# shellcheck source=/dev/null
source "${PROJECT_ROOT}/config/hpc_env.sh"
hpc_setup_step_04c_assembly_cleanup

export THREADS="${SLURM_CPUS_PER_TASK:-8}"
bash "${PROJECT_ROOT}/scripts/04c_clean_assemblies.sh"
