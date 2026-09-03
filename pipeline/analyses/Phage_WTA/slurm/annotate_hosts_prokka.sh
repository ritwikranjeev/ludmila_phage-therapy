#!/bin/bash
#SBATCH --job-name=HOST_ANN
#SBATCH --account=katrine_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/host_annotate_%x_%A_%a.out
#SBATCH --error=logs/host_annotate_%x_%A_%a.err

set -euo pipefail

if [[ "${DRY_RUN:-0}" == "1" ]]; then
  export PHAGE_WTA_SKIP_CONDA=1
fi

: "${SLURM_ARRAY_TASK_ID:?Submit with: sbatch --array=1-M slurm/annotate_hosts_prokka.sh (M = number of data rows in config/host_genomes.tsv)}"

# Slurm copies the batch script to the spool directory; paths from $0 break. Submit from project root.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  cd "${SLURM_SUBMIT_DIR}"
else
  cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
# shellcheck source=/dev/null
source "${PWD}/config/env_setup.sh"

cd "${PHAGE_WTA_PROJECT_ROOT}"
mkdir -p logs

TSV="${PHAGE_WTA_PROJECT_ROOT}/config/host_genomes.tsv"
# shellcheck source=/dev/null
source "${PHAGE_WTA_PROJECT_ROOT}/slurm/_array_tsv_common.sh"
array_tsv_validate_range "${TSV}" "${SLURM_ARRAY_TASK_ID}" "config/host_genomes.tsv" || exit 1

DATA_LINE=$((SLURM_ARRAY_TASK_ID + 1))
IFS=$'\t' read -r ID FASTA < <(sed -n "${DATA_LINE}p" "${TSV}")

if [[ -n "${FASTA:-}" && "${FASTA}" != /* ]]; then
  FASTA="${PHAGE_WTA_PROJECT_ROOT}/${FASTA}"
fi

if [[ -z "${ID:-}" || -z "${FASTA:-}" ]]; then
  echo "ERROR: No host genome row for array task ${SLURM_ARRAY_TASK_ID} (line ${DATA_LINE} of ${TSV})." >&2
  exit 1
fi

OUTDIR="${PHAGE_WTA_PROJECT_ROOT}/data/host_prokka/${ID}"
mkdir -p "${OUTDIR}"

hostname
date
echo "Prokka host: id=${ID} fasta=${FASTA} outdir=${OUTDIR}"

# Prokka with --outdir and --prefix writes PREFIX.faa, PREFIX.gff, PREFIX.tsv, etc.
# Expected proteins: "${OUTDIR}/${ID}.faa" (default naming).

# DRY_RUN=1: echo the command only. Export before sbatch, e.g.:
#   export DRY_RUN=1 && sbatch --array=1-1 slurm/annotate_hosts_prokka.sh
if [[ "${DRY_RUN:-0}" == "1" ]]; then
  echo "[DRY_RUN] Skipping prokka; would run:"
  echo prokka \
    --kingdom Bacteria \
    --cpus "${SLURM_CPUS_PER_TASK:-8}" \
    --outdir "${OUTDIR}" \
    --prefix "${ID}" \
    "${FASTA}"
  echo "[DRY_RUN] Expected output proteins: ${OUTDIR}/${ID}.faa"
  if [[ ! -f "${FASTA}" ]]; then
    echo "[DRY_RUN] NOTE: FASTA not found at ${FASTA} (real run would fail until path is valid)." >&2
  fi
  exit 0
fi

if [[ ! -f "${FASTA}" ]]; then
  echo "ERROR: FASTA not found: ${FASTA}" >&2
  exit 1
fi

prokka \
  --kingdom Bacteria \
  --cpus "${SLURM_CPUS_PER_TASK:-8}" \
  --outdir "${OUTDIR}" \
  --prefix "${ID}" \
  "${FASTA}"

echo "Done: ${OUTDIR}"
