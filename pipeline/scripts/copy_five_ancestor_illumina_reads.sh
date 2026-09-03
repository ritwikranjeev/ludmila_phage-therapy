#!/usr/bin/env bash
set -euo pipefail

# Copy five Illumina PE ancestor libraries from 3-res-mut raw_fastqs into
# data/raw/<short_name>/ (same layout as arbank sync).
#
# Edit SOURCE_DIR if your cluster path differs.
#
# On HPC3:
#   bash scripts/copy_five_ancestor_illumina_reads.sh

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SOURCE_DIR="${SOURCE_DIR:-/pub/ritwik1/manuscripts/ludmila-characterization/3-res-mut/reads/raw/raw_fastqs}"
DEST_ROOT="${PROJECT_ROOT}/data/raw"

die() { echo "ERROR: $*" >&2; exit 1; }
[[ -d "${SOURCE_DIR}" ]] || die "SOURCE_DIR not found: ${SOURCE_DIR}"

# dest_folder -> source basename prefix (files: ${prefix}_R1.fastq.gz and _R2)
declare -A MAP=(
  [PedCF]="S30.PedCF.Anc"
  [USA300]="S6.USA300.Anc"
  [OCPT04SA01]="S18.OCPT04SA01.Anc"
  [OCPT02SA01]="S12.OCPT02SA01.Anc"
  [SA412]="S24.SA412.Anc"
)

for dest in "${!MAP[@]}"; do
  pfx="${MAP[$dest]}"
  r1="${SOURCE_DIR}/${pfx}_R1.fastq.gz"
  r2="${SOURCE_DIR}/${pfx}_R2.fastq.gz"
  [[ -f "${r1}" && -f "${r2}" ]] || die "Missing ${r1} or ${r2}"
  out="${DEST_ROOT}/${dest}"
  mkdir -p "${out}"
  rsync -av "${r1}" "${r2}" "${out}/"
  echo "OK ${dest} <- ${pfx}_R{1,2}.fastq.gz"
done

echo "Done. Folders under ${DEST_ROOT}: PedCF USA300 OCPT04SA01 OCPT02SA01 SA412"
