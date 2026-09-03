#!/usr/bin/env bash
set -euo pipefail

# Run QUAST across all assemblies.
# Submit with: sbatch slurm/04_quast.slurm

# Environment: slurm/*.slurm sources config/hpc_env.sh

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ASSEMBLY_ROOT="${PROJECT_ROOT}/outputs/assemblies"
OUT_DIR="${PROJECT_ROOT}/outputs/quast"
THREADS="${THREADS:-8}"

die() { echo "ERROR: $*" >&2; exit 1; }
mkdir -p "${OUT_DIR}"
[[ -d "${ASSEMBLY_ROOT}" ]] || die "Missing ${ASSEMBLY_ROOT}"

mapfile -t assemblies < <(ls "${ASSEMBLY_ROOT}"/*/contigs.fasta 2>/dev/null || true)
[[ ${#assemblies[@]} -gt 0 ]] || die "No assembly files found under ${ASSEMBLY_ROOT}"

quast.py \
  -o "${OUT_DIR}" \
  -t "${THREADS}" \
  "${assemblies[@]}"

echo "QUAST completed: ${OUT_DIR}"
