#!/usr/bin/env bash
set -euo pipefail

# Assembly by platform:
# - illumina: SPAdes (or Unicycler short mode)
# - long: Flye (PacBio by default; not ONT)
# - hybrid: Unicycler hybrid mode
# Submit with: sbatch slurm/03_assemble.slurm
#
# Long-read Flye: PacBio preset (override with FLYE_PACBIO_MODE before sbatch).
# BOR libraries failed with --pacbio-hifi (no overlaps; high read divergence in logs) but work as raw:
#   FLYE_PACBIO_MODE=raw   -> --pacbio-raw  (default)
#   FLYE_PACBIO_MODE=hifi  -> --pacbio-hifi (use only if reads are true HiFi / CCS)

# Environment: slurm/*.slurm sources config/hpc_env.sh

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SAMPLES_TSV="${PROJECT_ROOT}/config/samples.tsv"
RAW_ROOT="${PROJECT_ROOT}/data/raw"
QC_FASTP="${PROJECT_ROOT}/outputs/qc/fastp"
OUT_ROOT="${PROJECT_ROOT}/outputs/assemblies"
THREADS="${THREADS:-16}"
FLYE_PACBIO_MODE="${FLYE_PACBIO_MODE:-raw}"

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*"; }

# Flye's conda-installed launcher often has #!/dfs6/.../python (dead on HPC). Run via env python.
flye_run() {
  [[ -n "${CONDA_PREFIX:-}" ]] || die "CONDA_PREFIX unset (source config/hpc_env.sh + hpc_setup_step_03_assemble)"
  local py="${CONDA_PREFIX}/bin/python3.11"
  [[ -x "${py}" ]] || py="${CONDA_PREFIX}/bin/python"
  [[ -x "${py}" ]] || die "No python in ${CONDA_PREFIX}/bin"
  local launcher="${CONDA_PREFIX}/bin/flye"
  [[ -f "${launcher}" ]] || die "Missing ${launcher}"
  "${py}" "${launcher}" "$@"
}

mkdir -p "${OUT_ROOT}"
[[ -f "${SAMPLES_TSV}" ]] || die "Missing ${SAMPLES_TSV}"
mapfile -t ROWS < <(tail -n +2 "${SAMPLES_TSV}")
[[ ${#ROWS[@]} -gt 0 ]] || die "No samples in ${SAMPLES_TSV}"

start_idx=0
end_idx=$((${#ROWS[@]} - 1))
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  start_idx="${SLURM_ARRAY_TASK_ID}"
  end_idx="${SLURM_ARRAY_TASK_ID}"
fi

for i in $(seq "${start_idx}" "${end_idx}"); do
  [[ "${i}" -lt "${#ROWS[@]}" ]] || die "Array index ${i} out of range"
  IFS=$'\t' read -r sample platform r1 r2 longread <<<"${ROWS[$i]}"

  sample_out="${OUT_ROOT}/${sample}"
  mkdir -p "${sample_out}"

  # Prefer trimmed reads when available.
  r1_in="${QC_FASTP}/${sample}/${sample}_R1.trimmed.fastq.gz"
  r2_in="${QC_FASTP}/${sample}/${sample}_R2.trimmed.fastq.gz"
  [[ -s "${r1_in}" ]] || r1_in="${RAW_ROOT}/${sample}/$(basename "${r1}")"
  [[ -s "${r2_in}" ]] || r2_in="${RAW_ROOT}/${sample}/$(basename "${r2}")"
  long_in="${RAW_ROOT}/${sample}/$(basename "${longread}")"

  log "Assembling ${sample} (${platform})"
  case "${platform}" in
    illumina)
      [[ -s "${r1_in}" && -s "${r2_in}" ]] || die "Missing Illumina reads for ${sample}"
      spades.py \
        -1 "${r1_in}" \
        -2 "${r2_in}" \
        -o "${sample_out}/spades" \
        --threads "${THREADS}"
      cp -f "${sample_out}/spades/contigs.fasta" "${sample_out}/contigs.fasta"
      ;;
    long)
      [[ -s "${long_in}" ]] || die "Missing long reads for ${sample}"
      case "${FLYE_PACBIO_MODE}" in
        raw)  flye_lr=(--pacbio-raw "${long_in}") ;;
        hifi) flye_lr=(--pacbio-hifi "${long_in}") ;;
        *) die "FLYE_PACBIO_MODE must be raw or hifi, got: ${FLYE_PACBIO_MODE}" ;;
      esac
      flye_run \
        "${flye_lr[@]}" \
        --out-dir "${sample_out}/flye" \
        --threads "${THREADS}"
      cp -f "${sample_out}/flye/assembly.fasta" "${sample_out}/contigs.fasta"
      ;;
    hybrid)
      [[ -s "${r1_in}" && -s "${r2_in}" && -s "${long_in}" ]] || die "Missing hybrid inputs for ${sample}"
      unicycler \
        -1 "${r1_in}" \
        -2 "${r2_in}" \
        -l "${long_in}" \
        -o "${sample_out}/unicycler" \
        --threads "${THREADS}"
      cp -f "${sample_out}/unicycler/assembly.fasta" "${sample_out}/contigs.fasta"
      ;;
    *)
      die "Unknown platform '${platform}' for ${sample}"
      ;;
  esac

  [[ -s "${sample_out}/contigs.fasta" ]] || die "Assembly output missing for ${sample}"
  log "Assembly completed for ${sample}"
done
