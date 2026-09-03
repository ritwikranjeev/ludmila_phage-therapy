#!/usr/bin/env bash
set -euo pipefail

# Prophage detection (default: PhiSpy). Override: PROPHAGE_TOOL=other
# Submit: sbatch slurm/11_prophage_scan.slurm
#
# PhiSpy 5 expects GenBank and runs SeqIO as genbank; use Prokka output
# outputs/annotation/<sample>/<sample>.gbk. LOCUS lines from Prokka are rewritten
# on the fly (see scripts/fix_prokka_gbk_locus.py) so Biopython accepts them.

# Environment: slurm/*.slurm sources config/hpc_env.sh

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SAMPLES_TSV="${PROJECT_ROOT}/config/samples.tsv"
ANNOT_ROOT="${PROJECT_ROOT}/outputs/annotation"
OUT_ROOT="${PROJECT_ROOT}/outputs/prophages"
PROPHAGE_TOOL="${PROPHAGE_TOOL:-phispy}"  # phispy|other
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-2}}"
FIX_GBK="${PROJECT_ROOT}/scripts/fix_prokka_gbk_locus.py"

mapfile -t ROWS < <(tail -n +2 "${SAMPLES_TSV}")
start=0; end=$((${#ROWS[@]} - 1))
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then start="${SLURM_ARRAY_TASK_ID}"; end="${SLURM_ARRAY_TASK_ID}"; fi
mkdir -p "${OUT_ROOT}"

for i in $(seq "${start}" "${end}"); do
  IFS=$'\t' read -r sample _ _ _ _ <<<"${ROWS[$i]}"
  gbk="${ANNOT_ROOT}/${sample}/${sample}.gbk"
  [[ -s "${gbk}" ]] || {
    echo "Missing Prokka GenBank for ${sample}: ${gbk} (run 05_prokka first)" >&2
    exit 1
  }
  work="${OUT_ROOT}/${sample}"
  mkdir -p "${work}"
  fixed="${work}/.input_fixed.gbk"
  python3 "${FIX_GBK}" "${gbk}" > "${fixed}"

  echo "PhiSpy ${sample}: ${gbk}" >&2
  if [[ "${PROPHAGE_TOOL}" == "phispy" ]]; then
    # output_choice 1 = prophage_coordinates.tsv only (default 3 also writes GBK;
    # SeqIO then errors on Prokka records missing molecule_type in annotations).
    if command -v phispy >/dev/null 2>&1; then
      phispy "${fixed}" -o "${work}" --threads "${THREADS}" --quiet --output_choice 1
    elif command -v PhiSpy.py >/dev/null 2>&1; then
      PhiSpy.py "${fixed}" -o "${work}" --threads "${THREADS}" --quiet --output_choice 1
    else
      echo "Neither phispy nor PhiSpy.py on PATH after conda_activate" >&2
      exit 1
    fi
  else
    echo "Please add '${PROPHAGE_TOOL}' invocation for ${sample}" > "${work}/TOOL_PLACEHOLDER.txt"
  fi
  rm -f "${fixed}"
done
