#!/usr/bin/env bash
set -euo pipefail

# Run pyseer on Panaroo gene presence/absence.
# Submit: sbatch slurm/14_pyseer_gene.slurm

# Environment: slurm/*.slurm sources config/hpc_env.sh

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PHENO="${PROJECT_ROOT}/config/phenotypes.tsv"
GENE_MATRIX="${PROJECT_ROOT}/outputs/pangenome/gene_presence_absence.Rtab"
OUT_DIR="${PROJECT_ROOT}/outputs/gwas/pyseer_gene"
# shellcheck source=scripts/ensure_population_dist.sh
source "${PROJECT_ROOT}/scripts/ensure_population_dist.sh"

mkdir -p "${OUT_DIR}"
ensure_population_dist

for trait in binary infection_class strong_infect; do
  pyseer \
    --phenotypes "${PHENO}" \
    --phenotype-column "${trait}" \
    --pres "${GENE_MATRIX}" \
    --distances "${PYSEER_STRUCTURE_DIST}" \
    > "${OUT_DIR}/pyseer_${trait}.txt"
done
