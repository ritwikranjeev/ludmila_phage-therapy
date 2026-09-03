#!/usr/bin/env bash
set -euo pipefail

# Read-level contamination diagnostics (USA300 mapping, Staph vs E. coli panel, contig depth+GC).
# Targets samples listed in config/diagnose_contamination_samples.txt (default: BOR-10, BOR-08, BOR-15).
#
# Submit: sbatch slurm/05a_diagnose_contamination.slurm
# Uses original outputs/assemblies/<sample>/contigs.fasta (not cleaned).

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LIST="${DIAGNOSE_SAMPLE_LIST:-${PROJECT_ROOT}/config/diagnose_contamination_samples.txt}"
SAMPLES_TSV="${PROJECT_ROOT}/config/samples.tsv"
ASSEMBLY_ROOT="${PROJECT_ROOT}/outputs/assemblies"
OUT_DIR="${PROJECT_ROOT}/outputs/diagnose_contamination"
USA300_REF="${PROJECT_ROOT}/data/reference/CP000046.1.fasta"
PANEL_REF="${PROJECT_ROOT}/data/reference/panel_Staph_Ecoli.fasta"
THREADS="${THREADS:-8}"

die() { echo "ERROR: $*" >&2; exit 1; }
[[ -f "${LIST}" ]] || die "Missing ${LIST}"
[[ -f "${SAMPLES_TSV}" ]] || die "Missing ${SAMPLES_TSV}"

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  mapfile -t SAMPLES < <(grep -vE '^[[:space:]]*($|#)' "${LIST}")
  [[ ${#SAMPLES[@]} -gt 0 ]] || die "No samples in ${LIST}"
  if (( SLURM_ARRAY_TASK_ID >= ${#SAMPLES[@]} )); then
    echo "Task ${SLURM_ARRAY_TASK_ID} out of range — exit 0"
    exit 0
  fi
  sample="${SAMPLES[${SLURM_ARRAY_TASK_ID}]}"
  asm="${ASSEMBLY_ROOT}/${sample}/contigs.fasta"
  [[ -s "${asm}" ]] || die "Missing ${asm}"
  python3 "${PROJECT_ROOT}/scripts/diagnose_contamination_one.py" \
    --sample "${sample}" \
    --samples-tsv "${SAMPLES_TSV}" \
    --assembly "${asm}" \
    --usa300-ref "${USA300_REF}" \
    --panel-ref "${PANEL_REF}" \
    --out-dir "${OUT_DIR}" \
    --threads "${THREADS}"
else
  while IFS= read -r line || [[ -n "${line}" ]]; do
    [[ -z "${line}" || "${line}" =~ ^[[:space:]]*# ]] && continue
    sample="${line//[[:space:]]/}"
    [[ -n "${sample}" ]] || continue
    asm="${ASSEMBLY_ROOT}/${sample}/contigs.fasta"
    [[ -s "${asm}" ]] || die "Missing ${asm}"
    python3 "${PROJECT_ROOT}/scripts/diagnose_contamination_one.py" \
      --sample "${sample}" \
      --samples-tsv "${SAMPLES_TSV}" \
      --assembly "${asm}" \
      --usa300-ref "${USA300_REF}" \
      --panel-ref "${PANEL_REF}" \
      --out-dir "${OUT_DIR}" \
      --threads "${THREADS}"
  done < "${LIST}"

  agg="${OUT_DIR}/all_samples_diagnose_summary.tsv"
  {
    echo -e "sample\tkey\tvalue"
    for d in "${OUT_DIR}"/*/; do
      [[ -f "${d}/diagnose_summary.txt" ]] || continue
      s="$(basename "${d}")"
      while IFS=$'\t' read -r k v; do
        [[ -z "${k}" ]] && continue
        echo -e "${s}\t${k}\t${v}"
      done < "${d}/diagnose_summary.txt"
    done
  } > "${agg}"
  echo "Wrote ${agg}"
fi
