#!/usr/bin/env bash
set -euo pipefail

# Run read QC by platform:
# - Illumina/hybrid: FastQC + fastp
# - Long/hybrid: NanoPlot
# Submit with: sbatch slurm/02_qc_reads.slurm

# Environment: slurm/*.slurm sources config/hpc_env.sh

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SAMPLES_TSV="${PROJECT_ROOT}/config/samples.tsv"
RAW_ROOT="${PROJECT_ROOT}/data/raw"
OUT="${PROJECT_ROOT}/outputs"
QC_FASTQC="${OUT}/qc/fastqc"
QC_FASTP="${OUT}/qc/fastp"
QC_NANOPLOT="${OUT}/qc/long_read"
THREADS="${THREADS:-8}"

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*"; }

mkdir -p "${QC_FASTQC}" "${QC_FASTP}" "${QC_NANOPLOT}"
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
  sample_raw="${RAW_ROOT}/${sample}"
  sample_fastqc="${QC_FASTQC}/${sample}"
  sample_fastp="${QC_FASTP}/${sample}"
  sample_nanoplot="${QC_NANOPLOT}/${sample}"
  mkdir -p "${sample_fastqc}" "${sample_fastp}" "${sample_nanoplot}"

  log "QC for ${sample} (${platform})"

  if [[ "${platform}" == "illumina" || "${platform}" == "hybrid" ]]; then
    r1_local="${sample_raw}/$(basename "${r1}")"
    r2_local="${sample_raw}/$(basename "${r2}")"
    [[ -s "${r1_local}" && -s "${r2_local}" ]] || die "Illumina reads missing for ${sample}"

    fastqc -t "${THREADS}" -o "${sample_fastqc}" "${r1_local}" "${r2_local}"

    fastp \
      -i "${r1_local}" \
      -I "${r2_local}" \
      -o "${sample_fastp}/${sample}_R1.trimmed.fastq.gz" \
      -O "${sample_fastp}/${sample}_R2.trimmed.fastq.gz" \
      -h "${sample_fastp}/${sample}.fastp.html" \
      -j "${sample_fastp}/${sample}.fastp.json" \
      -w "${THREADS}"
  fi

  if [[ "${platform}" == "long" || "${platform}" == "hybrid" ]]; then
    lr_local="${sample_raw}/$(basename "${longread}")"
    [[ -s "${lr_local}" ]] || die "Long read file missing for ${sample}"
    # FastQC works on any FASTQ (Illumina or long-read); cheap baseline QC.
    fastqc -t "${THREADS}" -o "${sample_fastqc}" "${lr_local}"
    # NanoPlot (conda env) — optional if binary missing on PATH
    if command -v NanoPlot >/dev/null 2>&1; then
      NanoPlot \
        --fastq "${lr_local}" \
        --threads "${THREADS}" \
        --outdir "${sample_nanoplot}" \
        --prefix "${sample}_"
    else
      log "WARN: NanoPlot not on PATH; long-read QC is FastQC-only for ${sample}"
    fi
  fi

  log "QC completed for ${sample}"
done
