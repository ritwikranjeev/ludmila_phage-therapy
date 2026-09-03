#!/bin/bash
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH -t 48:00:00
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -J snakemake_mrsa_ont
#SBATCH -o /pub/ritwik1/manuscripts/OCPT-Ludmila_Characterization/analyses/MRSA_ONT_comparative/logs/snakemake_%j.out
#SBATCH -e /pub/ritwik1/manuscripts/OCPT-Ludmila_Characterization/analyses/MRSA_ONT_comparative/logs/snakemake_%j.err
#SBATCH --export=ALL

# Snakemake driver: submits per-rule jobs to SLURM. Check logs/snakemake_<jobid>.out

set -euo pipefail
ANALYSIS_DIR="/pub/ritwik1/manuscripts/OCPT-Ludmila_Characterization/analyses/MRSA_ONT_comparative"
cd "$ANALYSIS_DIR"
mkdir -p logs

set +e
# Reuse Figure5 conda + snakemake env on HPC3
source "/pub/ritwik1/manuscripts/OCPT-Ludmila_Characterization/analyses/Figure5_Phage_resistant_Mutant_Analysis/scripts/activate_conda_snakemake.sh"
set -e

export XDG_CACHE_HOME="$ANALYSIS_DIR/.snakemake/xdg_cache"
mkdir -p "$XDG_CACHE_HOME"

CONDA_PREFIX_DIR="$ANALYSIS_DIR/.snakemake/conda"
mkdir -p "$CONDA_PREFIX_DIR"

echo "Running MRSA ONT comparative Snakemake from $(pwd) at $(date)"
snakemake --profile profile/slurm --use-conda \
  --conda-prefix "$CONDA_PREFIX_DIR" \
  --rerun-incomplete --latency-wait 180

echo "Snakemake finished at $(date)"
