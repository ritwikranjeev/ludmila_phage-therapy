#!/usr/bin/env bash
set -euo pipefail

# Run defense annotation (default: DefenseFinder). Override: DEFENSE_TOOL=padloc
# Submit: sbatch slurm/08_defense_annotation.slurm
#
# DefenseFinder nucleotide path uses pyhmmer; some envs hit seq.name.decode() on str (pyhmmer update).
# Default for defensefinder: use Prokka proteins outputs/annotation/<sample>/<sample>.faa when present
# (skips prodigal/pyhmmer DNA reader). DEFENSE_INPUT_MODE=nucleotide forces contigs.fasta instead.
# PADLOC always uses nucleotide contigs (cleaned preferred like Prokka).
# DEFENSE_USE_ORIGINAL_ASSEMBLIES=1 skips cleaned assemblies for nucleotide input.

# Environment: slurm/*.slurm sources config/hpc_env.sh

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SAMPLES_TSV="${PROJECT_ROOT}/config/samples.tsv"
ASSEMBLY_ROOT="${PROJECT_ROOT}/outputs/assemblies"
ASSEMBLY_CLEANED="${PROJECT_ROOT}/outputs/assemblies_cleaned"
ANNOT_ROOT="${PROJECT_ROOT}/outputs/annotation"
OUT_ROOT="${PROJECT_ROOT}/outputs/defense/raw"
DEFENSE_TOOL="${DEFENSE_TOOL:-defensefinder}"  # defensefinder|padloc
DEFENSE_INPUT_MODE="${DEFENSE_INPUT_MODE:-auto}" # auto|nucleotide|protein (defensefinder only)
THREADS="${THREADS:-8}"

mapfile -t ROWS < <(tail -n +2 "${SAMPLES_TSV}")
start=0; end=$((${#ROWS[@]} - 1))
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then start="${SLURM_ARRAY_TASK_ID}"; end="${SLURM_ARRAY_TASK_ID}"; fi
mkdir -p "${OUT_ROOT}"

for i in $(seq "${start}" "${end}"); do
  IFS=$'\t' read -r sample _ _ _ _ <<<"${ROWS[$i]}"
  cleaned="${ASSEMBLY_CLEANED}/${sample}/contigs.fasta"
  orig="${ASSEMBLY_ROOT}/${sample}/contigs.fasta"
  faa="${ANNOT_ROOT}/${sample}/${sample}.faa"
  fasta=""
  if [[ "${DEFENSE_TOOL}" == "defensefinder" ]]; then
    if [[ "${DEFENSE_INPUT_MODE}" == "protein" ]]; then
      [[ -s "${faa}" ]] || { echo "DEFENSE_INPUT_MODE=protein but missing ${faa}" >&2; exit 1; }
      fasta="${faa}"
    elif [[ "${DEFENSE_INPUT_MODE}" == "nucleotide" ]]; then
      if [[ -z "${DEFENSE_USE_ORIGINAL_ASSEMBLIES:-}" ]] && [[ -s "${cleaned}" ]]; then
        fasta="${cleaned}"
      elif [[ -s "${orig}" ]]; then
        fasta="${orig}"
      else
        echo "Missing assembly for ${sample}" >&2
        exit 1
      fi
    else
      # auto: prefer Prokka FAA (avoids DefenseFinder DNA/pyhmmer reader bugs)
      if [[ -s "${faa}" ]]; then
        fasta="${faa}"
      elif [[ -z "${DEFENSE_USE_ORIGINAL_ASSEMBLIES:-}" ]] && [[ -s "${cleaned}" ]]; then
        fasta="${cleaned}"
      elif [[ -s "${orig}" ]]; then
        fasta="${orig}"
      else
        echo "Missing input for ${sample} (no ${faa} or assemblies)" >&2
        exit 1
      fi
    fi
  else
    # padloc: nucleotide only
    if [[ -z "${DEFENSE_USE_ORIGINAL_ASSEMBLIES:-}" ]] && [[ -s "${cleaned}" ]]; then
      fasta="${cleaned}"
    elif [[ -s "${orig}" ]]; then
      fasta="${orig}"
    else
      echo "Missing assembly for ${sample} (tried ${cleaned} and ${orig})" >&2
      exit 1
    fi
  fi
  echo "Defense ${sample}: ${fasta}" >&2
  out="${OUT_ROOT}/${sample}"
  mkdir -p "${out}"
  if [[ "${DEFENSE_TOOL}" == "defensefinder" ]]; then
    defense-finder run "${fasta}" -o "${out}" --workers "${THREADS}"
  elif [[ "${DEFENSE_TOOL}" == "padloc" ]]; then
    padloc --fna "${fasta}" --outdir "${out}" --cpu "${THREADS}"
  else
    echo "Unsupported DEFENSE_TOOL=${DEFENSE_TOOL}" >&2
    exit 1
  fi
done
