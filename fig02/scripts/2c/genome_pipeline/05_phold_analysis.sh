#!/bin/bash
# Run Phold on Pharokka GenBank output (second step of Pharokka → Phold pipeline).
# Input: results/annotation/pharokka/pharokka.gbk (run pharokka_annotation.sh first)
# Output: results/annotation/phold/

#SBATCH --job-name=phage_phold
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
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

if [[ ! -f "$PHAROKKA_GBK" ]]; then
  echo "ERROR: Pharokka GenBank not found. Run pharokka_annotation.sh first. Expected $PHAROKKA_GBK" >&2
  exit 1
fi

mkdir -p "$PHOLD_DIR"

"$PHOLD_PYTHON" "$PHOLD_BIN" predict -i "$PHAROKKA_GBK" -o "$PHOLD_DIR" -f -t 16

echo "Done. Phold output: $PHOLD_DIR"
