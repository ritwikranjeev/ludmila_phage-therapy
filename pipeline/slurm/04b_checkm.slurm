#!/usr/bin/env bash
# Submit: sbatch slurm/04b_checkm.slurm
# CheckM taxonomy_wf (genus Staphylococcus) — completeness / contamination for all assemblies.

#SBATCH -J checkm
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 08:00:00
#SBATCH -o outputs/logs/%x_%A.out
#SBATCH -e outputs/logs/%x_%A.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR}"
# shellcheck source=/dev/null
source "${PROJECT_ROOT}/config/hpc_env.sh"
hpc_setup_step_04b_checkm

bash "${PROJECT_ROOT}/scripts/04b_checkm.sh"
