#!/bin/bash
# Annotate phage assembly with Pharokka (first step of Pharokka → Phold pipeline).
# Input: results/assembly/{PHAGE_NAME}.fasta (run assembly.sh first)
# Output: results/annotation/pharokka/ (includes pharokka.gbk for phold)

#SBATCH --job-name=phage_pharokka
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

if [[ ! -f "$ASSEMBLY_FASTA" ]]; then
  echo "ERROR: Assembly not found. Run assembly.sh first. Expected $ASSEMBLY_FASTA" >&2
  exit 1
fi

if [[ -z "$PHAROKKA_DB" ]]; then
  echo "ERROR: Set PHAROKKA_DB in config/config.sh to the Pharokka database directory." >&2
  exit 1
fi

mkdir -p "$PHAROKKA_DIR"

"$PHAROKKA_PYTHON" "$PHAROKKA_BIN" -i "$ASSEMBLY_FASTA" -o "$PHAROKKA_DIR" -f -d "$PHAROKKA_DB" -t 16

echo "Done. Pharokka output: $PHAROKKA_DIR (use pharokka.gbk for phold_analysis.sh)"
