#!/usr/bin/env bash
set -euo pipefail

# Build manifest files after Prokka:
# - config/gff_list.txt
# - config/sample_to_paths.tsv

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ANNOT_ROOT="${PROJECT_ROOT}/outputs/annotation"
ASSEMBLY_ROOT="${PROJECT_ROOT}/outputs/assemblies"
ASSEMBLY_CLEANED="${PROJECT_ROOT}/outputs/assemblies_cleaned"
GFF_LIST="${PROJECT_ROOT}/config/gff_list.txt"
MAP_TSV="${PROJECT_ROOT}/config/sample_to_paths.tsv"

echo -n > "${GFF_LIST}"
echo -e "sample\tgff\tcontigs" > "${MAP_TSV}"

for gff in "${ANNOT_ROOT}"/*/*.gff; do
  [[ -e "${gff}" ]] || continue
  sample="$(basename "${gff%.gff}")"
  cleaned="${ASSEMBLY_CLEANED}/${sample}/contigs.fasta"
  if [[ -s "${cleaned}" ]]; then
    contigs="${cleaned}"
  else
    contigs="${ASSEMBLY_ROOT}/${sample}/contigs.fasta"
  fi
  echo "${gff}" >> "${GFF_LIST}"
  echo -e "${sample}\t${gff}\t${contigs}" >> "${MAP_TSV}"
done

echo "Wrote ${GFF_LIST} and ${MAP_TSV}"
