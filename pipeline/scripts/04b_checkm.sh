#!/usr/bin/env bash
set -euo pipefail

# CheckM taxonomy workflow (completeness / contamination) for isolate assemblies.
# Uses genus-level markers for Staphylococcus (appropriate for this cohort).
# Submit with: sbatch slurm/04b_checkm.slurm
#
# Note: lineage_wf (genome tree + pplacer) is not used here — pplacer often fails
# on login nodes / low memory; taxonomy_wf avoids that and still reports QA metrics.

# Environment: slurm/*.slurm sources config/hpc_env.sh

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ASSEMBLY_ROOT="${PROJECT_ROOT}/outputs/assemblies"
CHECKM_ROOT="${PROJECT_ROOT}/outputs/checkm"
BINS_DIR="${CHECKM_ROOT}/bins"
OUT_DIR="${CHECKM_ROOT}/taxonomy_staphylococcus"
THREADS="${THREADS:-16}"
RANK="${CHECKM_TAXON_RANK:-genus}"
TAXON="${CHECKM_TAXON:-Staphylococcus}"

die() { echo "ERROR: $*" >&2; exit 1; }

[[ -d "${ASSEMBLY_ROOT}" ]] || die "Missing ${ASSEMBLY_ROOT}"
mapfile -t fastas < <(ls "${ASSEMBLY_ROOT}"/*/contigs.fasta 2>/dev/null || true)
[[ ${#fastas[@]} -gt 0 ]] || die "No contigs.fasta under ${ASSEMBLY_ROOT}"

mkdir -p "${BINS_DIR}"
rm -rf "${OUT_DIR}"

for f in "${fastas[@]}"; do
  sample="$(basename "$(dirname "${f}")")"
  ln -sf "$(realpath "${f}")" "${BINS_DIR}/${sample}.fasta"
done

echo "CheckM bins: ${#fastas[@]} genomes in ${BINS_DIR}"
checkm taxonomy_wf "${RANK}" "${TAXON}" -t "${THREADS}" -x fasta "${BINS_DIR}" "${OUT_DIR}"

echo "CheckM completed: ${OUT_DIR}"
echo "Main table: ${OUT_DIR}/storage/bin_stats_ext.tsv"
