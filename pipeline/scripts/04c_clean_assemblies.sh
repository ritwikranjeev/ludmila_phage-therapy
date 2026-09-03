#!/usr/bin/env bash
set -euo pipefail

# Clean high-contamination assemblies: read depth vs contigs + map contigs to USA300 chr (CP000046.1).
# Outputs: outputs/assemblies_cleaned/<sample>/contigs.fasta and per-sample removed-contig logs.
#
# Submit (parallel, 4 at a time): sbatch slurm/04c_clean_assemblies.slurm
# Run all listed samples in one shell (login/debug): bash scripts/04c_clean_assemblies.sh

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SAMPLES_FILE="${ASSEMBLY_CLEANUP_LIST:-${PROJECT_ROOT}/config/assembly_cleanup_samples.txt}"
SAMPLES_TSV="${PROJECT_ROOT}/config/samples.tsv"
ASSEMBLY_ROOT="${PROJECT_ROOT}/outputs/assemblies"
OUT_ROOT="${PROJECT_ROOT}/outputs/assemblies_cleaned"
REF_FASTA="${REFERENCE_SAUREUS_CHR:-${PROJECT_ROOT}/data/reference/CP000046.1.fasta}"
THREADS="${THREADS:-8}"

die() { echo "ERROR: $*" >&2; exit 1; }
[[ -f "${SAMPLES_FILE}" ]] || die "Missing ${SAMPLES_FILE}"
[[ -f "${SAMPLES_TSV}" ]] || die "Missing ${SAMPLES_TSV}"
[[ -d "${ASSEMBLY_ROOT}" ]] || die "Missing ${ASSEMBLY_ROOT}"

run_one() {
  local sample="$1"
  local contigs="${ASSEMBLY_ROOT}/${sample}/contigs.fasta"
  [[ -s "${contigs}" ]] || die "Missing assembly for ${sample}: ${contigs}"

  local outdir="${OUT_ROOT}/${sample}"
  mkdir -p "${outdir}"
  python3 "${PROJECT_ROOT}/scripts/assembly_cleanup_one.py" \
    --sample "${sample}" \
    --contigs "${contigs}" \
    --samples-tsv "${SAMPLES_TSV}" \
    --ref "${REF_FASTA}" \
    --out-fasta "${outdir}/contigs.fasta" \
    --out-removed "${outdir}/removed_contigs.tsv" \
    --threads "${THREADS}"
}

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  mapfile -t LIST < <(grep -vE '^[[:space:]]*($|#)' "${SAMPLES_FILE}")
  [[ ${#LIST[@]} -gt 0 ]] || die "No samples in ${SAMPLES_FILE}"
  if (( SLURM_ARRAY_TASK_ID >= ${#LIST[@]} )); then
    echo "Array task ${SLURM_ARRAY_TASK_ID} >= ${#LIST[@]} — nothing to do."
    exit 0
  fi
  sample="${LIST[${SLURM_ARRAY_TASK_ID}]}"
  echo "Array task ${SLURM_ARRAY_TASK_ID}: ${sample}"
  run_one "${sample}"
else
  while IFS= read -r line || [[ -n "${line}" ]]; do
    [[ -z "${line}" || "${line}" =~ ^[[:space:]]*# ]] && continue
    sample="${line//[[:space:]]/}"
    [[ -n "${sample}" ]] || continue
    run_one "${sample}"
  done < "${SAMPLES_FILE}"
  echo "All cleanup jobs finished → ${OUT_ROOT}"
fi
