#!/usr/bin/env bash
set -euo pipefail

# Annotate assemblies with Prokka.
# Submit with:
#   sbatch slurm/05_prokka.slurm          — all samples in config/samples.tsv (array 0–82)
#   sbatch slurm/05_prokka_cleaned_only.slurm — only samples listed in config/assembly_cleanup_samples.txt (11 jobs)
#
# Assembly FASTA: prefers outputs/assemblies_cleaned/<sample>/contigs.fasta when present
# (post 04c cleanup); otherwise outputs/assemblies/<sample>/contigs.fasta.
# Set PROKKA_USE_ORIGINAL_ASSEMBLIES=1 to skip cleaned and always use outputs/assemblies.
#
# To annotate an arbitrary list: PROKKA_SAMPLE_LIST=/path/to/file.txt sbatch --array=0-N ...

# Environment: slurm/*.slurm sources config/hpc_env.sh

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SAMPLES_TSV="${PROJECT_ROOT}/config/samples.tsv"
ASSEMBLY_ROOT="${PROJECT_ROOT}/outputs/assemblies"
ASSEMBLY_CLEANED="${PROJECT_ROOT}/outputs/assemblies_cleaned"
ANNOT_ROOT="${PROJECT_ROOT}/outputs/annotation"
THREADS="${THREADS:-8}"
# Default list for “cleaned-cohort only” is the same 11 strains as assembly cleanup.
PROKKA_SAMPLE_LIST_DEFAULT="${PROJECT_ROOT}/config/assembly_cleanup_samples.txt"

die(){ echo "ERROR: $*" >&2; exit 1; }

run_prokka_for_sample() {
  local sample="$1"
  local cleaned="${ASSEMBLY_CLEANED}/${sample}/contigs.fasta"
  local orig="${ASSEMBLY_ROOT}/${sample}/contigs.fasta"
  local fasta=""
  if [[ -z "${PROKKA_USE_ORIGINAL_ASSEMBLIES:-}" ]] && [[ -s "${cleaned}" ]]; then
    fasta="${cleaned}"
  elif [[ -s "${orig}" ]]; then
    fasta="${orig}"
  else
    die "Missing assembly for ${sample} (tried ${cleaned} and ${orig})"
  fi
  local out="${ANNOT_ROOT}/${sample}"
  echo "Prokka ${sample}: ${fasta}" >&2
  prokka \
    --force \
    --outdir "${out}" \
    --prefix "${sample}" \
    --locustag "${sample//-/}" \
    --genus Staphylococcus \
    --species aureus \
    --cpus "${THREADS}" \
    "${fasta}"
}

mkdir -p "${ANNOT_ROOT}"

LIST_FILE="${PROKKA_SAMPLE_LIST:-}"
if [[ -z "${LIST_FILE}" && -n "${PROKKA_CLEANED_ONLY:-}" ]]; then
  LIST_FILE="${PROKKA_SAMPLE_LIST_DEFAULT}"
fi

if [[ -n "${LIST_FILE}" ]]; then
  [[ -f "${LIST_FILE}" ]] || die "PROKKA_SAMPLE_LIST not found: ${LIST_FILE}"
  mapfile -t NAMES < <(
    grep -vE '^[[:space:]]*($|#)' "${LIST_FILE}" | sed 's/[[:space:]]//g' | grep -v '^$'
  )
  [[ ${#NAMES[@]} -gt 0 ]] || die "No sample names in ${LIST_FILE}"

  if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    if (( SLURM_ARRAY_TASK_ID >= ${#NAMES[@]} )); then
      echo "Array task ${SLURM_ARRAY_TASK_ID} >= ${#NAMES[@]} — nothing to do."
      exit 0
    fi
    run_prokka_for_sample "${NAMES[$SLURM_ARRAY_TASK_ID]}"
  else
    for sample in "${NAMES[@]}"; do
      run_prokka_for_sample "${sample}"
    done
  fi
  exit 0
fi

mapfile -t ROWS < <(tail -n +2 "${SAMPLES_TSV}")

start=0; end=$((${#ROWS[@]} - 1))
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then start="${SLURM_ARRAY_TASK_ID}"; end="${SLURM_ARRAY_TASK_ID}"; fi

for i in $(seq "${start}" "${end}"); do
  IFS=$'\t' read -r sample _ _ _ _ <<<"${ROWS[$i]}"
  run_prokka_for_sample "${sample}"
done
