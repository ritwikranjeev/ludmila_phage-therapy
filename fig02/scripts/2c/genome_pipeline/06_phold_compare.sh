#!/bin/bash
# Optional: Compare Phold and Pharokka annotations (which genes Phold added, agreement).
# Input: Pharokka .gbk and Phold predictions dir (run pharokka_annotation.sh and phold_analysis.sh first)
# Output: results/annotation/phold/compare/

#SBATCH --job-name=phage_phold_compare
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
# compare invokes foldseek; Slurm PATH often lacks pholdENV/bin (foldseek lives next to phold).
export PATH="$(dirname "$PHOLD_PYTHON"):${PATH}"

if [[ ! -f "$PHAROKKA_GBK" ]]; then
  echo "ERROR: Pharokka GenBank not found. Run pharokka_annotation.sh first. Expected $PHAROKKA_GBK" >&2
  exit 1
fi

mkdir -p "$PHOLD_COMPARE_DIR"

"$PHOLD_PYTHON" "$PHOLD_BIN" compare -i "$PHAROKKA_GBK" --predictions_dir "$PHOLD_DIR" -o "$PHOLD_COMPARE_DIR" -f -t 8

echo "Done. Compare output: $PHOLD_COMPARE_DIR"
