#!/usr/bin/env bash
set -euo pipefail

# QUAST on cleaned assemblies only (outputs/assemblies_cleaned/<sample>/contigs.fasta).
# Stages uniquely named symlinks so QUAST labels are per-sample (not all "contigs").
# Submit with: sbatch slurm/04d_quast_cleaned.slurm

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CLEANED_ROOT="${PROJECT_ROOT}/outputs/assemblies_cleaned"
OUT_DIR="${PROJECT_ROOT}/outputs/quast_cleaned"
STAGING="${OUT_DIR}/staging"
THREADS="${THREADS:-8}"

die() { echo "ERROR: $*" >&2; exit 1; }
[[ -d "${CLEANED_ROOT}" ]] || die "Missing ${CLEANED_ROOT} — run 04c cleanup first"

mkdir -p "${STAGING}"
rm -f "${STAGING}"/*.fasta

mapfile -t samples < <(
  shopt -s nullglob
  for d in "${CLEANED_ROOT}"/*/; do
    [[ -d "${d}" ]] && basename "${d}"
  done | sort
)
[[ ${#samples[@]} -gt 0 ]] || die "No sample directories under ${CLEANED_ROOT}"

inputs=()
for sample in "${samples[@]}"; do
  f="${CLEANED_ROOT}/${sample}/contigs.fasta"
  [[ -s "${f}" ]] || continue
  ln -sf "$(realpath "${f}")" "${STAGING}/${sample}.fasta"
  inputs+=("${STAGING}/${sample}.fasta")
done
[[ ${#inputs[@]} -gt 0 ]] || die "No contigs.fasta found under ${CLEANED_ROOT}"

quast.py \
  -o "${OUT_DIR}" \
  -t "${THREADS}" \
  "${inputs[@]}"

echo "QUAST (cleaned) completed: ${OUT_DIR}"
