#!/bin/bash
# Clean raw Illumina reads with fastp (same parameters as Figure 5 pipeline).
# Input: data/raw_fastqs/{PHAGE_NAME}_R1.fastq.gz, _R2.fastq.gz
# Output: results/reads/clean/, results/qc/fastp/

#SBATCH --job-name=phage_clean_reads
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --error=%x.%A.err
#SBATCH --out=%x.%A.out
#SBATCH --mail-type=fail,end
#SBATCH --mail-user=ritwik1@uci.edu

set -e
# When run via sbatch, SLURM runs from spool dir; use submit dir as project root. Otherwise use script location.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  cd "$SLURM_SUBMIT_DIR"
else
  cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/.."
fi
source config/config.sh

module load fastp/0.23.2 2>/dev/null || true

R1="$RAW_READS_DIR/${PHAGE_NAME}_R1.fastq.gz"
R2="$RAW_READS_DIR/${PHAGE_NAME}_R2.fastq.gz"
if [[ ! -f "$R1" || ! -f "$R2" ]]; then
  echo "ERROR: Raw reads not found. Expected $R1 and $R2" >&2
  exit 1
fi

mkdir -p "$READS_CLEAN" "$QC_DIR/fastp"

fastp \
  -i "$R1" -I "$R2" \
  -o "$READS_CLEAN/${PHAGE_NAME}_clean_R1.fq.gz" -O "$READS_CLEAN/${PHAGE_NAME}_clean_R2.fq.gz" \
  --detect_adapter_for_pe \
  --cut_front --cut_tail \
  --cut_window_size 4 --cut_mean_quality 20 \
  --length_required 60 \
  --correction \
  --thread 16 \
  -h "$QC_DIR/fastp/fastp.html" -j "$QC_DIR/fastp/fastp.json"

echo "Done. Clean reads: $READS_CLEAN/${PHAGE_NAME}_clean_R1.fq.gz and _R2.fq.gz"
