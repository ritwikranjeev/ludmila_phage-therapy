#!/usr/bin/env bash
set -euo pipefail

# CheckM taxonomy_wf on cleaned assemblies (same genus markers as 04b_checkm).
# Submit with: sbatch slurm/04e_checkm_cleaned.slurm

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CLEANED_ROOT="${PROJECT_ROOT}/outputs/assemblies_cleaned"
CHECKM_ROOT="${PROJECT_ROOT}/outputs/checkm_cleaned"
BINS_DIR="${CHECKM_ROOT}/bins"
OUT_DIR="${CHECKM_ROOT}/taxonomy_staphylococcus"
THREADS="${THREADS:-16}"
RANK="${CHECKM_TAXON_RANK:-genus}"
TAXON="${CHECKM_TAXON:-Staphylococcus}"

die() { echo "ERROR: $*" >&2; exit 1; }
[[ -d "${CLEANED_ROOT}" ]] || die "Missing ${CLEANED_ROOT} — run 04c cleanup first"

mapfile -t samples < <(
  shopt -s nullglob
  for d in "${CLEANED_ROOT}"/*/; do
    [[ -d "${d}" ]] && basename "${d}"
  done | sort
)
[[ ${#samples[@]} -gt 0 ]] || die "No sample directories under ${CLEANED_ROOT}"

mkdir -p "${BINS_DIR}"
rm -rf "${OUT_DIR}"

count=0
for sample in "${samples[@]}"; do
  f="${CLEANED_ROOT}/${sample}/contigs.fasta"
  [[ -s "${f}" ]] || continue
  ln -sf "$(realpath "${f}")" "${BINS_DIR}/${sample}.fasta"
  count=$((count + 1))
done
[[ "${count}" -gt 0 ]] || die "No contigs.fasta under ${CLEANED_ROOT}"

echo "CheckM (cleaned): ${count} genomes in ${BINS_DIR}"
checkm taxonomy_wf "${RANK}" "${TAXON}" -t "${THREADS}" -x fasta "${BINS_DIR}" "${OUT_DIR}"

echo "CheckM (cleaned) completed: ${OUT_DIR}"
echo "Main table: ${OUT_DIR}/storage/bin_stats_ext.tsv"
