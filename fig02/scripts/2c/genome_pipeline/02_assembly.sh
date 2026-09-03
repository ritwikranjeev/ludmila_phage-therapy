#!/bin/bash
# Assemble phage genome from cleaned reads using SPAdes --metaviral.
# Input: results/reads/clean/{PHAGE_NAME}_clean_R1.fq.gz, _R2.fq.gz (run clean_reads.sh first)
# Output: results/assembly/{PHAGE_NAME}.fasta, results/qc/seqkit, optional results/qc/quast

#SBATCH --job-name=phage_assembly
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
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

module load SPAdes/3.15.4 2>/dev/null || true
# Seqkit: use conda base or module if available
if [[ -n "${CONDA_INIT:-}" && -f "${CONDA_INIT}" ]]; then
  source "$CONDA_INIT"
  conda activate base 2>/dev/null || true
fi

R1="$READS_CLEAN/${PHAGE_NAME}_clean_R1.fq.gz"
R2="$READS_CLEAN/${PHAGE_NAME}_clean_R2.fq.gz"
if [[ ! -f "$R1" || ! -f "$R2" ]]; then
  echo "ERROR: Clean reads not found. Run clean_reads.sh first. Expected $R1 and $R2" >&2
  exit 1
fi

mkdir -p "$ASSEMBLY_DIR" "$QC_DIR/seqkit" "$QC_DIR/quast"

# SPAdes with --metaviral
spades.py \
  -1 "$R1" -2 "$R2" \
  --metaviral -t 16 -m 120 \
  -o "$ASSEMBLY_DIR"

# Use scaffolds.fasta or contigs.fasta (metaviral may keep output in K*/ or use different names)
SPADES_FASTA=""
for f in "$ASSEMBLY_DIR/scaffolds.fasta" "$ASSEMBLY_DIR/contigs.fasta"; do
  if [[ -f "$f" ]]; then
    SPADES_FASTA="$f"
    break
  fi
done
if [[ -z "$SPADES_FASTA" ]]; then
  for kdir in "$ASSEMBLY_DIR"/K*/; do
    if [[ -f "${kdir}scaffolds.fasta" ]]; then
      SPADES_FASTA="${kdir}scaffolds.fasta"
      break
    fi
    if [[ -f "${kdir}final_contigs.fasta" ]]; then
      SPADES_FASTA="${kdir}final_contigs.fasta"
      break
    fi
  done
fi
if [[ -z "$SPADES_FASTA" || ! -f "$SPADES_FASTA" ]]; then
  echo "ERROR: No assembly FASTA found in $ASSEMBLY_DIR (tried scaffolds.fasta, contigs.fasta, K*/scaffolds.fasta, K*/final_contigs.fasta)." >&2
  exit 1
fi

# Filter contigs >= 500 bp
seqkit seq -m 500 "$SPADES_FASTA" > "$ASSEMBLY_FASTA"

# Seqkit stats
seqkit stats "$SPADES_FASTA" "$ASSEMBLY_FASTA" > "$QC_DIR/seqkit/seqkit.stats.txt"

# Optional: QUAST (set QUAST_PYTHON in config to path to quast.py to enable)
if [[ -n "${QUAST_PYTHON:-}" && -f "$QUAST_PYTHON" ]]; then
  python3 "$QUAST_PYTHON" "$ASSEMBLY_FASTA" -o "$QC_DIR/quast" -t 16
fi

echo "Done. Assembly: $ASSEMBLY_FASTA"
