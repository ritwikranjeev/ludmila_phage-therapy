#!/usr/bin/env bash
set -euo pipefail

# Run Panaroo on all Prokka GFFs.
# Submit with: sbatch slurm/06_panaroo.slurm

# Environment: slurm/*.slurm sources config/hpc_env.sh

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GFF_LIST="${PROJECT_ROOT}/config/gff_list.txt"
OUT_DIR="${PROJECT_ROOT}/outputs/pangenome"
THREADS="${THREADS:-16}"

[[ -s "${GFF_LIST}" ]] || { echo "Missing ${GFF_LIST}. Run scripts/make_annotation_manifests.sh first."; exit 1; }
mkdir -p "${OUT_DIR}"

# Panaroo puts tempfile.mkdtemp() inside -o; on NFS (errno 70 / joblib "Communication error on send")
# is common under heavy parallelism. Use node-local scratch on Slurm (TMPDIR may be unset — fall back to /tmp).
PANAROO_OUT="${OUT_DIR}"
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  _local="${SLURM_TMPDIR:-/tmp}"
  PANAROO_OUT="${_local}/panaroo_${SLURM_JOB_ID}"
  rm -rf "${PANAROO_OUT}"
  mkdir -p "${PANAROO_OUT}"
  echo "Panaroo writing to node-local: ${PANAROO_OUT}"
fi

panaroo \
  -i $(tr '\n' ' ' < "${GFF_LIST}") \
  -o "${PANAROO_OUT}" \
  --clean-mode strict \
  -t "${THREADS}" \
  -a core

if [[ "${PANAROO_OUT}" != "${OUT_DIR}" ]]; then
  mkdir -p "${OUT_DIR}"
  rsync -a --delete "${PANAROO_OUT}/" "${OUT_DIR}/"
  rm -rf "${PANAROO_OUT}"
  echo "Synced pangenome results to ${OUT_DIR}"
fi
