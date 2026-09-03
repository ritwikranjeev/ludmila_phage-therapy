#!/bin/bash
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH -t 24:00:00
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -J breseq_S16
#SBATCH -o logs/breseq_S16_%j.out
#SBATCH -e logs/breseq_S16_%j.err

# One-off breseq for S16.OCPT04SA01.4.6 (avoids Snakemake worker I/O on home).
# After this succeeds, run: sbatch run_snakemake.sh  (it will skip breseq and do gdtools + tables + plot).

set -e
ANALYSIS_DIR="/pub/ritwik1/manuscripts/OCPT-Ludmila_Characterization/analyses/Figure5_Phage_resistant_Mutant_Analysis"
cd "$ANALYSIS_DIR"
mkdir -p logs

# Use project dir for cache so we don't hit I/O error on /data/homezvol1/$USER
export XDG_CACHE_HOME="$ANALYSIS_DIR/.snakemake/xdg_cache"
mkdir -p "$XDG_CACHE_HOME"

# Activate conda then the breseq env (same as pipeline)
CONDA_SH=""
for d in miniconda3_snakemake miniconda3 miniforge3; do
  if [[ -f "$HOME/$d/etc/profile.d/conda.sh" ]]; then
    CONDA_SH="$HOME/$d/etc/profile.d/conda.sh"
    break
  fi
done
if [[ -z "$CONDA_SH" ]]; then
  echo "ERROR: No conda profile.d/conda.sh found." >&2
  exit 1
fi
source "$CONDA_SH"
conda activate "$ANALYSIS_DIR/envs/conda/ae711351d26df9fc6a32b273dbacf033_"

REF_GBK="$ANALYSIS_DIR/results/references/S18.OCPT04SA01.Anc/S18.OCPT04SA01.Anc.gbk"
OUT_DIR="$ANALYSIS_DIR/results/breseq/S16.OCPT04SA01.4.6"
R1="$ANALYSIS_DIR/results/qc/cleaned/S16.OCPT04SA01.4.6_R1.fq.gz"
R2="$ANALYSIS_DIR/results/qc/cleaned/S16.OCPT04SA01.4.6_R2.fq.gz"

for f in "$REF_GBK" "$R1" "$R2"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: Missing input $f" >&2
    exit 1
  fi
done

mkdir -p "$OUT_DIR"
echo "Running breseq for S16 at $(date)"
breseq -r "$REF_GBK" -o "$OUT_DIR" -j 16 --name S16.OCPT04SA01.4.6 "$R1" "$R2"
echo "Breseq S16 finished at $(date)"

if [[ -f "$OUT_DIR/output/output.gd" ]]; then
  echo "SUCCESS: $OUT_DIR/output/output.gd created. Run: sbatch run_snakemake.sh"
else
  echo "WARNING: output.gd not found. Check logs/breseq_S16_$SLURM_JOB_ID.err"
  exit 1
fi
