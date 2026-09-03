#!/usr/bin/env bash
set -euo pipefail

# Copy reads listed in config/samples.tsv into data/raw/<sample>/.
# Submit with: sbatch slurm/01_fetch_reads.slurm
# Local run:   bash scripts/01_fetch_reads.sh
# Note for interactive debug:
#   srun -A katrine_lab --pty bash -i

# Environment: slurm/*.slurm sources config/hpc_env.sh

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SAMPLES_TSV="${PROJECT_ROOT}/config/samples.tsv"
DEST_ROOT="${PROJECT_ROOT}/data/raw"
COPY_CMD="${COPY_CMD:-rsync}"  # options: rsync|cp

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*"; }

[[ -f "${SAMPLES_TSV}" ]] || die "Missing ${SAMPLES_TSV}"
mkdir -p "${DEST_ROOT}"

mapfile -t ROWS < <(tail -n +2 "${SAMPLES_TSV}")
[[ ${#ROWS[@]} -gt 0 ]] || die "No samples in ${SAMPLES_TSV}"

start_idx=0
end_idx=$((${#ROWS[@]} - 1))
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  start_idx="${SLURM_ARRAY_TASK_ID}"
  end_idx="${SLURM_ARRAY_TASK_ID}"
fi

copy_one() {
  local src="$1" dst_dir="$2"
  [[ "${src}" == "NA" || -z "${src}" ]] && return 0
  [[ -f "${src}" ]] || die "Source file not found: ${src}"
  if [[ "${COPY_CMD}" == "rsync" ]]; then
    rsync -avh --progress "${src}" "${dst_dir}/"
  elif [[ "${COPY_CMD}" == "cp" ]]; then
    cp -v "${src}" "${dst_dir}/"
  else
    die "Unsupported COPY_CMD=${COPY_CMD}; use rsync or cp"
  fi
}

for i in $(seq "${start_idx}" "${end_idx}"); do
  [[ "${i}" -lt "${#ROWS[@]}" ]] || die "Array index ${i} out of range"
  IFS=$'\t' read -r sample platform r1 r2 longread <<<"${ROWS[$i]}"
  [[ -n "${sample}" ]] || die "Empty sample value at row $((i + 2))"
  sample_dir="${DEST_ROOT}/${sample}"
  mkdir -p "${sample_dir}"
  log "Copying reads for ${sample} (${platform})"

  copy_one "${r1}" "${sample_dir}"
  copy_one "${r2}" "${sample_dir}"
  copy_one "${longread}" "${sample_dir}"

  for maybe_file in "${r1}" "${r2}" "${longread}"; do
    [[ "${maybe_file}" == "NA" || -z "${maybe_file}" ]] && continue
    copied="${sample_dir}/$(basename "${maybe_file}")"
    [[ -s "${copied}" ]] || die "Copied file missing/empty: ${copied}"
  done

  log "Completed ${sample}"
done

log "Fetch step completed."
