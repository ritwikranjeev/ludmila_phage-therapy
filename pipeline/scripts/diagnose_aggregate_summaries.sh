#!/usr/bin/env bash
# After array 05a finishes, merge per-sample diagnose_summary.txt into one TSV.
# Usage: bash scripts/diagnose_aggregate_summaries.sh

set -euo pipefail
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${PROJECT_ROOT}/outputs/diagnose_contamination"
agg="${OUT_DIR}/all_samples_diagnose_summary.tsv"
{
  echo -e "sample\tkey\tvalue"
  for f in "${OUT_DIR}"/*/diagnose_summary.txt; do
    [[ -f "${f}" ]] || continue
    s="$(basename "$(dirname "${f}")")"
    while IFS=$'\t' read -r k v; do
      [[ -z "${k}" ]] && continue
      echo -e "${s}\t${k}\t${v}"
    done < "${f}"
  done
} > "${agg}"
echo "Wrote ${agg}"
