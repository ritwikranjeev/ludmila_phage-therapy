#!/usr/bin/env bash
set -euo pipefail

# Copy per-strain read folders from the arbank SRA staging tree into this project:
#   data/raw/<sample>/*.fastq(.gz)
#
# Run on HPC3 (from project root or any cwd — resolves paths from this script):
#   bash scripts/copy_reads_from_arbank.sh
#
# Large copy in background + log:
#   mkdir -p logs
#   nohup bash scripts/copy_reads_from_arbank.sh >> logs/copy_reads_from_arbank.log 2>&1 &
#   tail -f logs/copy_reads_from_arbank.log
#
# Optional: dry run
#   DRY_RUN=1 bash scripts/copy_reads_from_arbank.sh
#
# Override source if your reads move:
#   SOURCE_ROOT=/path/to/data bash scripts/copy_reads_from_arbank.sh

# ---- Editable defaults (HPC3 / Ludmila arbank layout) ----
SOURCE_ROOT="${SOURCE_ROOT:-/pub/ritwik1/manuscripts/ludmila-characterization/1-host-range-phylogeny/arbank/info/sra/data}"
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST_ROOT="${DEST_ROOT:-${PROJECT_ROOT}/data/raw}"

log() { echo "[$(date -Is)] $*"; }

if [[ ! -d "${SOURCE_ROOT}" ]]; then
  echo "ERROR: SOURCE_ROOT does not exist: ${SOURCE_ROOT}" >&2
  exit 1
fi

mkdir -p "${DEST_ROOT}"

RSYNC_FLAGS=(-a)
if [[ "${DRY_RUN:-0}" == "1" ]]; then
  RSYNC_FLAGS+=(-n -v)
  log "DRY RUN (no files copied)"
fi

n=0
shopt -s nullglob
for src_dir in "${SOURCE_ROOT}"/*/; do
  [[ -d "${src_dir}" ]] || continue
  sample="$(basename "${src_dir}")"
  # Skip obvious non-sample dirs if someone adds them next to strain folders.
  case "${sample}" in
    toolkit|lost+found) continue ;;
  esac

  dest="${DEST_ROOT}/${sample}"
  mkdir -p "${dest}"
  log "Sync ${sample}  (${src_dir} -> ${dest}/)"
  rsync "${RSYNC_FLAGS[@]}" "${src_dir}" "${dest}/"
  n=$((n + 1))
done

log "Finished: ${n} sample folders synced under ${DEST_ROOT}"
