#!/usr/bin/env bash
set -euo pipefail

# pyseer on curated host WTA / receptor gene presence (Prokka panel).
# Submit: sbatch slurm/21_pyseer_host_receptors.slurm
# Prereqs: Prokka (05), phenotypes (13), IQ-TREE (07), scripts/19_extract_host_receptor_genes.py

# Environment: slurm/*.slurm sources config/hpc_env.sh

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PHENO="${PROJECT_ROOT}/config/phenotypes.tsv"
MATRIX="${PROJECT_ROOT}/outputs/receptor_candidates/host_receptor_presence.Rtab"
OUT_DIR="${PROJECT_ROOT}/outputs/gwas/pyseer_host_receptors"
# shellcheck source=scripts/ensure_population_dist.sh
source "${PROJECT_ROOT}/scripts/ensure_population_dist.sh"

mkdir -p "${OUT_DIR}"
[[ -s "${MATRIX}" ]] || {
  echo "Missing ${MATRIX}. Run: python3 scripts/19_extract_host_receptor_genes.py" >&2
  exit 1
}

ensure_population_dist

for trait in binary infection_class strong_infect; do
  pyseer \
    --phenotypes "${PHENO}" \
    --phenotype-column "${trait}" \
    --pres "${MATRIX}" \
    --distances "${PYSEER_STRUCTURE_DIST}" \
    > "${OUT_DIR}/pyseer_host_receptors_${trait}.txt"
done
