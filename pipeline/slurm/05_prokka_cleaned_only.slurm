#!/usr/bin/env bash
# Submit: sbatch slurm/05_prokka_cleaned_only.slurm
# Prokka only for strains in config/assembly_cleanup_samples.txt (11 samples; same list as 04c cleanup).
# Array 0–10, max 20 concurrent (only 11 tasks will run).

#SBATCH -J prokka-cln
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH --array=0-10%20
#SBATCH -o outputs/logs/%x_%A_%a.out
#SBATCH -e outputs/logs/%x_%A_%a.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR}"
# shellcheck source=/dev/null
source "${PROJECT_ROOT}/config/hpc_env.sh"
hpc_setup_step_05_prokka

export PROKKA_CLEANED_ONLY=1
export PROKKA_SAMPLE_LIST="${PROKKA_SAMPLE_LIST:-${PROJECT_ROOT}/config/assembly_cleanup_samples.txt}"
export THREADS="${SLURM_CPUS_PER_TASK:-8}"

bash "${PROJECT_ROOT}/scripts/05_prokka.sh"
