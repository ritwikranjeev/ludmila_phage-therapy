#!/bin/bash
# Source from other Slurm drivers:  source "${PHAGE_WTA_PROJECT_ROOT}/slurm/_array_tsv_common.sh"
# Validates SLURM_ARRAY_TASK_ID against data row count in a TSV (line 1 = header).

array_tsv_count_data_rows() {
  local tsv="$1"
  if [[ ! -f "$tsv" ]]; then
    echo "ERROR: TSV not found: ${tsv}" >&2
    return 1
  fi
  local n
  n=$(wc -l < "$tsv" | tr -d ' ' | head -c 20)
  if [[ -z "$n" || "$n" -lt 2 ]]; then
    echo "ERROR: ${tsv} must have a header line and at least one data row." >&2
    return 1
  fi
  echo $((n - 1))
}

# Usage: array_tsv_validate_range /path/to.tsv "$SLURM_ARRAY_TASK_ID" "label"
array_tsv_validate_range() {
  local tsv="$1"
  local task="$2"
  local label="${3:-TSV}"
  local nd
  nd=$(array_tsv_count_data_rows "$tsv") || return 1
  if [[ "$task" -lt 1 ]] || [[ "$task" -gt "$nd" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID=${task} out of range for ${label} (valid 1-${nd} per data rows in ${tsv})." >&2
    return 1
  fi
  return 0
}
