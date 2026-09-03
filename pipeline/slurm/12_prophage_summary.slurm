#!/usr/bin/env bash
# Submit after step 11 completes: sbatch --dependency=afterok:<11_jobid> slurm/12_prophage_summary.slurm
# Or: sbatch slurm/12_prophage_summary.slurm
#SBATCH -J prophage_summary
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH -o outputs/logs/%x_%A.out
#SBATCH -e outputs/logs/%x_%A.err
set -euo pipefail
PROJECT_ROOT="${SLURM_SUBMIT_DIR}"
# shellcheck source=/dev/null
source "${PROJECT_ROOT}/config/hpc_env.sh"
hpc_setup_step_12_prophage_summary

python3 "${PROJECT_ROOT}/scripts/12_summarize_prophages.py"
