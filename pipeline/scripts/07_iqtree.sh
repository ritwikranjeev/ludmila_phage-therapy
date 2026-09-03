#!/usr/bin/env bash
set -euo pipefail

# Build core-genome ML tree from Panaroo core alignment.
# Submit with: sbatch slurm/07_iqtree.slurm

# Environment: slurm/*.slurm sources config/hpc_env.sh

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ALIGN="${PROJECT_ROOT}/outputs/pangenome/core_gene_alignment.aln"
OUT_DIR="${PROJECT_ROOT}/outputs/trees"
THREADS="${THREADS:-16}"

mkdir -p "${OUT_DIR}"
[[ -s "${ALIGN}" ]] || { echo "Missing alignment: ${ALIGN}"; exit 1; }

# Conda iqtree2 package often installs the binary as `iqtree` (IQ-TREE 2).
IQTREE_CMD="$(command -v iqtree2 2>/dev/null || command -v iqtree 2>/dev/null || true)"
[[ -n "${IQTREE_CMD}" ]] || {
  echo "Neither iqtree2 nor iqtree on PATH after activating CONDA_ENV_IQTREE" >&2
  exit 1
}

"${IQTREE_CMD}" \
  -s "${ALIGN}" \
  -m GTR+G \
  -B 1000 \
  -alrt 1000 \
  -T "${THREADS}" \
  --prefix "${OUT_DIR}/core"
