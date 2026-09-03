#!/bin/bash
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH -t 48:00:00
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -J snakemake_fig5
#SBATCH -o logs/snakemake_%j.out
#SBATCH -e logs/snakemake_%j.err
#SBATCH --export=ALL

# Snakemake driver job: runs on one node and submits workflow jobs via SLURM.
# You can log out after sbatch; check progress with squeue -u $USER and logs/snakemake_<jobid>.out

set -e
ANALYSIS_DIR="/pub/ritwik1/manuscripts/OCPT-Ludmila_Characterization/analyses/Figure5_Phage_resistant_Mutant_Analysis"
cd "$ANALYSIS_DIR"
mkdir -p logs

# Activate conda (HPC3: snakemake env is often under pub housekeeping miniforge)
set +e
source "$ANALYSIS_DIR/scripts/activate_conda_snakemake.sh"
set -e
# Use project dir for Snakemake cache so compute nodes don't hit I/O errors on home (/data/homezvol1/$USER)
export XDG_CACHE_HOME="$ANALYSIS_DIR/.snakemake/xdg_cache"
mkdir -p "$XDG_CACHE_HOME"

echo "Running Snakemake from $(pwd) at $(date)"
snakemake -s workflow/Snakefile --profile profile/slurm --use-conda \
  --conda-prefix envs/conda --rerun-incomplete --rerun-triggers mtime --latency-wait 120 \
  --config project_dir="$ANALYSIS_DIR"
echo "Snakemake finished at $(date)"
