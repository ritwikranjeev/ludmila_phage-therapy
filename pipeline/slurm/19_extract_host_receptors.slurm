#!/usr/bin/env bash
# Run after Prokka + config/samples.tsv. Builds host WTA / receptor feature tables.
# Submit: sbatch slurm/19_extract_host_receptors.slurm

#SBATCH -J host_receptors
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
hpc_setup_step_19_host_receptors

python3 "${PROJECT_ROOT}/scripts/19_extract_host_receptor_genes.py"
