#!/bin/bash
#SBATCH --job-name=WTA_SCAN
#SBATCH --account=katrine_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/wta_scan_%x_%A_%a.out
#SBATCH --error=logs/wta_scan_%x_%A_%a.err

set -euo pipefail

# Stdlib-only driver; phage_wta_env may live on a filesystem not mounted on all compute nodes.
export PHAGE_WTA_SKIP_CONDA=1

: "${SLURM_ARRAY_TASK_ID:?Submit with: sbatch --array=1-M slurm/scan_hosts_for_wta_genes.sh (M = data rows in config/host_genomes.tsv)}"

# Slurm copies the batch script to the spool directory; paths from $0 break. Submit from project root.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  cd "${SLURM_SUBMIT_DIR}"
else
  cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
# shellcheck source=/dev/null
source "${PWD}/config/env_setup.sh"

cd "${PHAGE_WTA_PROJECT_ROOT}"
mkdir -p logs data/wta_features/by_host

TSV="${PHAGE_WTA_PROJECT_ROOT}/config/host_genomes.tsv"
# shellcheck source=/dev/null
source "${PHAGE_WTA_PROJECT_ROOT}/slurm/_array_tsv_common.sh"
array_tsv_validate_range "${TSV}" "${SLURM_ARRAY_TASK_ID}" "config/host_genomes.tsv" || exit 1

DATA_LINE=$((SLURM_ARRAY_TASK_ID + 1))
IFS=$'\t' read -r ID _FASTA < <(sed -n "${DATA_LINE}p" "${TSV}")

if [[ -z "${ID:-}" ]]; then
  echo "ERROR: No host row for array task ${SLURM_ARRAY_TASK_ID} (line ${DATA_LINE} of ${TSV})." >&2
  exit 1
fi

hostname
date
echo "WTA scan host_id=${ID}"

# DRY_RUN=1: echo only (export DRY_RUN=1 before sbatch).
if [[ "${DRY_RUN:-0}" == "1" ]]; then
  echo "[DRY_RUN] Would run:"
  echo "  python3 \"${PHAGE_WTA_PROJECT_ROOT}/bin/scan_hosts_for_wta_genes.py\" --host-id \"${ID}\""
  echo "[DRY_RUN] Expects host proteins: ${PHAGE_WTA_PROJECT_ROOT}/data/host_prokka/${ID}/${ID}.faa (or single *.faa)"
  echo "[DRY_RUN] HMM DB: ${PHAGE_WTA_PROJECT_ROOT}/data/wta_hmms/wta_combined.hmm (+ hmmpress indices)"
  exit 0
fi

python3 "${PHAGE_WTA_PROJECT_ROOT}/bin/scan_hosts_for_wta_genes.py" --host-id "${ID}"

echo "Done host ${ID}."
