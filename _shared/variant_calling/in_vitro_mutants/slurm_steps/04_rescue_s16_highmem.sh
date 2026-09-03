#!/bin/bash
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH -t 24:00:00
# Standard MaxMemPerCPU=6144MB: need enough CPUs so mem cap >= requested RAM.
# 30*6144 MiB ≈ 180 GiB cap; node ~191000 MiB — stay under with 180G request.
#SBATCH -c 30
#SBATCH --mem=180G
#SBATCH -J breseq_S16_rescue
#SBATCH -o logs/04_rescue_s16_%j.out
#SBATCH -e logs/04_rescue_s16_%j.err

# Dedicated rescue run for S16 (known high-memory sample).
# Uses low breseq parallelism (-j 1) and a fresh output directory.

set -euo pipefail
ANALYSIS_DIR="/pub/ritwik1/manuscripts/OCPT-Ludmila_Characterization/analyses/Figure5_Phage_resistant_Mutant_Analysis"
cd "$ANALYSIS_DIR"
mkdir -p logs

# Snakemake env lacks bowtie2/R; breseq needs the Snakemake-created breseq conda env on PATH.
BRESEQ_ENV="$ANALYSIS_DIR/envs/conda/ae711351d26df9fc6a32b273dbacf033_"
if [[ ! -x "$BRESEQ_ENV/bin/breseq" ]]; then
  echo "ERROR: breseq env missing: $BRESEQ_ENV/bin/breseq" >&2
  exit 1
fi
export PATH="$BRESEQ_ENV/bin:$PATH"
BRESEQ_BIN="$BRESEQ_ENV/bin/breseq"

REF_GBK="$ANALYSIS_DIR/results/references/S18.OCPT04SA01.Anc/S18.OCPT04SA01.Anc.gbk"
OUT_DIR="$ANALYSIS_DIR/results/breseq/S16.OCPT04SA01.4.6"
R1="$ANALYSIS_DIR/results/qc/cleaned/S16.OCPT04SA01.4.6_R1.fq.gz"
R2="$ANALYSIS_DIR/results/qc/cleaned/S16.OCPT04SA01.4.6_R2.fq.gz"

for f in "$REF_GBK" "$R1" "$R2"; do
  [[ -f "$f" ]] || { echo "ERROR: Missing input $f" >&2; exit 1; }
done

# Remove stale intermediate state so breseq cannot short-circuit on old files.
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

echo "Running S16 breseq rescue at $(date)"
echo "Node: $(hostname)"
"$BRESEQ_BIN" -r "$REF_GBK" -o "$OUT_DIR" -j 1 --name S16.OCPT04SA01.4.6 "$R1" "$R2"

echo "Finished at $(date)"
[[ -f "$OUT_DIR/output/output.gd" ]] || { echo "ERROR: output.gd missing" >&2; exit 1; }
echo "SUCCESS: $OUT_DIR/output/output.gd"
