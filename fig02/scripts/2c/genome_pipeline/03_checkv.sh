#!/bin/bash
# Assess phage genome quality with CheckV.
# Input: results/assembly/{PHAGE_NAME}.fasta (run assembly.sh first)
# Output: results/checkv/

#SBATCH --job-name=phage_checkv
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --error=%x.%A.err
#SBATCH --out=%x.%A.out
#SBATCH --mail-type=fail,end
#SBATCH --mail-user=ritwik1@uci.edu

set -e
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  cd "$SLURM_SUBMIT_DIR"
else
  cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/.."
fi
source config/config.sh

if [[ ! -f "$ASSEMBLY_FASTA" ]]; then
  echo "ERROR: Assembly not found. Run assembly.sh first. Expected $ASSEMBLY_FASTA" >&2
  exit 1
fi

# CheckV database: set CHECKVDB and/or pass -d
if [[ -n "$CHECKV_DB" ]]; then
  export CHECKVDB="$CHECKV_DB"
fi

mkdir -p "$CHECKV_DIR"

CHECKV_OPTS=(-t 16)
[[ -n "$CHECKV_DB" ]] && CHECKV_OPTS+=(-d "$CHECKV_DB")
"$CHECKV_PYTHON" "$CHECKV_BIN" end_to_end "$ASSEMBLY_FASTA" "$CHECKV_DIR" "${CHECKV_OPTS[@]}"
echo "Done. CheckV results: $CHECKV_DIR"
