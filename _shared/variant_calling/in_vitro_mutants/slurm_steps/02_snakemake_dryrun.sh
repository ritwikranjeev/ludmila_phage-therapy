#!/bin/bash
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH -t 02:00:00
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -J fig5_dryrun
#SBATCH -o logs/02_dryrun_%j.out
#SBATCH -e logs/02_dryrun_%j.err
#SBATCH --export=ALL
#
# Step 2: Snakemake --dry-run with SLURM profile (no compute-heavy rules executed).
# Confirms DAG and job plan. Run after step 1 succeeds.

set -euo pipefail
ANALYSIS_DIR="/pub/ritwik1/manuscripts/OCPT-Ludmila_Characterization/analyses/Figure5_Phage_resistant_Mutant_Analysis"
cd "$ANALYSIS_DIR"
mkdir -p logs

set +e
source "$ANALYSIS_DIR/scripts/activate_conda_snakemake.sh"
set -e

export XDG_CACHE_HOME="$ANALYSIS_DIR/.snakemake/xdg_cache"
mkdir -p "$XDG_CACHE_HOME"

echo "=== Figure5 step 2: snakemake --dry-run ==="
snakemake -s workflow/Snakefile --profile profile/slurm --dry-run --use-conda \
  --conda-prefix envs/conda --config project_dir="$ANALYSIS_DIR"
echo "=== dry-run finished ==="
