#!/bin/bash
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH -t 48:00:00
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -J fig5_snakemake
#SBATCH -o logs/03_full_pipeline_%j.out
#SBATCH -e logs/03_full_pipeline_%j.err
#SBATCH --export=ALL
#
# Step 3: Full Figure5a pipeline (Snakemake submits per-rule SLURM jobs).
# Driver walltime must cover the whole workflow (_refs, fastp, breseq, plots).
# Same logic as run_snakemake.sh in pipeline root; kept here for numbered workflow.

set -euo pipefail
ANALYSIS_DIR="/pub/ritwik1/manuscripts/OCPT-Ludmila_Characterization/analyses/Figure5_Phage_resistant_Mutant_Analysis"
cd "$ANALYSIS_DIR"
mkdir -p logs

set +e
source "$ANALYSIS_DIR/scripts/activate_conda_snakemake.sh"
set -e

export XDG_CACHE_HOME="$ANALYSIS_DIR/.snakemake/xdg_cache"
mkdir -p "$XDG_CACHE_HOME"

echo "=== Figure5 step 3: full Snakemake run ==="
echo "Started: $(date)"
snakemake -s workflow/Snakefile --profile profile/slurm --use-conda \
  --conda-prefix envs/conda --rerun-incomplete --rerun-triggers mtime --latency-wait 120 \
  --config project_dir="$ANALYSIS_DIR"
echo "Finished: $(date)"
