#!/usr/bin/env bash
# Sync this analysis folder to HPC3. By default excludes config/host_genomes.tsv so a laptop
# toy TSV does not overwrite the server list (re-run sync_host_genomes_from_host_prokka.py on HPC if needed).
#
# Usage (from anywhere):
#   bash bin/rsync_phage_wta_to_hpc.sh
#
# Env (optional):
#   PHAGE_WTA_HPC3_USER          default: ritwik1
#   PHAGE_WTA_HPC3_HOST          default: hpc3.rcic.uci.edu
#   PHAGE_WTA_HPC3_DEST          default: .../OCPT-Ludmila_Characterization/analyses/Phage_WTA/
#   PHAGE_WTA_RSYNC_INCLUDE_HOST_TSV=1   also sync config/host_genomes.tsv

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PHAGE_WTA_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

USER="${PHAGE_WTA_HPC3_USER:-ritwik1}"
HOST="${PHAGE_WTA_HPC3_HOST:-hpc3.rcic.uci.edu}"
DEST="${PHAGE_WTA_HPC3_DEST:-${USER}@${HOST}:/pub/ritwik1/manuscripts/OCPT-Ludmila_Characterization/analyses/Phage_WTA/}"

RSYNC_OPTS=(-avz --progress)
if [[ "${PHAGE_WTA_RSYNC_INCLUDE_HOST_TSV:-0}" != "1" ]]; then
  RSYNC_OPTS+=(--exclude=config/host_genomes.tsv)
fi
RSYNC_OPTS+=("${PHAGE_WTA_ROOT}/" "${DEST}")

rsync "${RSYNC_OPTS[@]}"
