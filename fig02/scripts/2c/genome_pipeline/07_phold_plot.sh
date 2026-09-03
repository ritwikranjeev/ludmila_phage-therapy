#!/bin/bash
# Optional: Circos plot of annotated genome from Phold compare output.
# Input: results/annotation/phold/compare/phold.gbk (run phold_compare.sh first)
# Output: results/annotation/phold/plot/

#SBATCH --job-name=phage_phold_plot
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

PHOLD_GBK="$PHOLD_COMPARE_DIR/phold.gbk"
if [[ ! -f "$PHOLD_GBK" ]]; then
  echo "ERROR: Phold compare GenBank not found. Run phold_compare.sh first. Expected $PHOLD_GBK" >&2
  exit 1
fi

mkdir -p "$PHOLD_PLOT_DIR"

"$PHOLD_PYTHON" "$PHOLD_BIN" plot -i "$PHOLD_GBK" -o "$PHOLD_PLOT_DIR" -t "$PHOLD_PLOT_TITLE" -f

echo "Done. Plot output: $PHOLD_PLOT_DIR"
