#!/bin/bash
#SBATCH --job-name=PHAGE_RBP
#SBATCH --account=katrine_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/phage_rbp_%x_%A_%a.out
#SBATCH --error=logs/phage_rbp_%x_%A_%a.err

set -euo pipefail

if [[ "${DRY_RUN:-0}" == "1" ]]; then
  export PHAGE_WTA_SKIP_CONDA=1
fi

: "${SLURM_ARRAY_TASK_ID:?Submit with: sbatch --array=1-N slurm/detect_phage_rbp.sh (N = number of data rows in config/phage_genomes.tsv)}"

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

TSV="${PHAGE_WTA_PROJECT_ROOT}/config/phage_genomes.tsv"
# shellcheck source=/dev/null
source "${PHAGE_WTA_PROJECT_ROOT}/slurm/_array_tsv_common.sh"
array_tsv_validate_range "${TSV}" "${SLURM_ARRAY_TASK_ID}" "config/phage_genomes.tsv" || exit 1

DATA_LINE=$((SLURM_ARRAY_TASK_ID + 1))
IFS=$'\t' read -r ID FASTA < <(sed -n "${DATA_LINE}p" "${TSV}")

if [[ -n "${FASTA:-}" && "${FASTA}" != /* ]]; then
  FASTA="${PHAGE_WTA_PROJECT_ROOT}/${FASTA}"
fi

if [[ -z "${ID:-}" || -z "${FASTA:-}" ]]; then
  echo "ERROR: No phage genome row for array task ${SLURM_ARRAY_TASK_ID} (line ${DATA_LINE} of ${TSV})." >&2
  exit 1
fi

# Prokka default: --prefix ID in --outdir produces ID.faa (matches annotate_phage_prokka.sh).
# If you changed --prefix, adjust this path or glob data/phage_prokka/${ID}/*.faa.
PROTEINS="${PHAGE_WTA_PROJECT_ROOT}/data/phage_prokka/${ID}/${ID}.faa"
OUTDIR="${PHAGE_WTA_PROJECT_ROOT}/data/phage_rbp/${ID}"
mkdir -p "${OUTDIR}"

hostname
date
echo "Phage RBP: id=${ID} proteins=${PROTEINS} outdir=${OUTDIR}"

# DRY_RUN=1: echo tool invocation only (export DRY_RUN=1 before sbatch).
if [[ "${DRY_RUN:-0}" == "1" ]]; then
  echo "[DRY_RUN] Would require Prokka proteins at: ${PROTEINS}"
  echo "[DRY_RUN] Would run (after TODO: confirm CLI) one of:"
  echo "  PhageRBPdetection \"${PROTEINS}\" \"${OUTDIR}\""
  echo "  phageRBPdetection.py \"${PROTEINS}\" \"${OUTDIR}\""
  exit 0
fi

if [[ ! -f "${PROTEINS}" ]]; then
  echo "ERROR: Expected Prokka protein FASTA not found: ${PROTEINS}" >&2
  echo "Hint: check data/phage_prokka/${ID}/*.faa and align naming with Prokka output." >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# TODO: Replace the command below with the exact PhageRBPdetection invocation
# from your install (`which`, README, -h/--help). Repo layouts differ (e.g. v2
# _standalone.py, v3 _inference.py, v4 PhageRBPdetect_v4_inference.py); you may
# need `python /path/to/script.py` plus model paths, threads, etc.
# TODO: Confirm input/output flags so all results are written under OUTDIR (cd,
#   --outdir, or copy step as required).
# The two-argument example is a placeholder only — adjust before production use.
# ---------------------------------------------------------------------------
if command -v PhageRBPdetection >/dev/null 2>&1; then
  PhageRBPdetection "${PROTEINS}" "${OUTDIR}"
elif command -v phageRBPdetection.py >/dev/null 2>&1; then
  phageRBPdetection.py "${PROTEINS}" "${OUTDIR}"
else
  echo "ERROR: Neither PhageRBPdetection nor phageRBPdetection.py found in PATH." >&2
  exit 1
fi

echo "Done: ${OUTDIR}"
