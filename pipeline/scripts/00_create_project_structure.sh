#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash scripts/00_create_project_structure.sh
#
# Note for interactive debugging on this cluster:
#   srun -A katrine_lab --pty bash -i

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${PROJECT_ROOT}/outputs"

mkdir -p \
  "${PROJECT_ROOT}/config" \
  "${PROJECT_ROOT}/data/raw" \
  "${PROJECT_ROOT}/data/phage" \
  "${PROJECT_ROOT}/data/phenotypes" \
  "${OUT}/qc/fastqc" \
  "${OUT}/qc/fastp" \
  "${OUT}/qc/long_read" \
  "${OUT}/assemblies" \
  "${OUT}/quast" \
  "${OUT}/annotation" \
  "${OUT}/pangenome" \
  "${OUT}/trees" \
  "${OUT}/defense/raw" \
  "${OUT}/prophages" \
  "${OUT}/gwas/pyseer_gene" \
  "${OUT}/gwas/pyseer_defense" \
  "${OUT}/gwas/pyseer_prophage" \
  "${OUT}/gwas/pyseer_host_receptors" \
  "${OUT}/receptor_candidates" \
  "${OUT}/logs" \
  "${PROJECT_ROOT}/scripts" \
  "${PROJECT_ROOT}/slurm" \
  "${PROJECT_ROOT}/notebooks" \
  "${PROJECT_ROOT}/envs" \
  "${PROJECT_ROOT}/logs"

echo "Project structure created under: ${PROJECT_ROOT}"
echo "Pipeline outputs (assemblies, annotation, GWAS, …): ${OUT}/"
