#!/bin/bash
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH -t 00:15:00
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -J fig5_linkreads
#SBATCH -o logs/00_link_reads_%j.out
#SBATCH -e logs/00_link_reads_%j.err
#
# Step 0 (HPC3): Symlink ancestral + mutant FASTQs into data/Figure5_.../sequencing/raw
# so samples.tsv paths exist. Run once before 01_validate_samples.sh.

set -euo pipefail
ANALYSIS_DIR="/pub/ritwik1/manuscripts/OCPT-Ludmila_Characterization/analyses/Figure5_Phage_resistant_Mutant_Analysis"
cd "$ANALYSIS_DIR"
mkdir -p logs

export PROJECT="/pub/ritwik1/manuscripts/OCPT-Ludmila_Characterization"
bash scripts/link_figure5_reads_hpc3.sh
echo "=== Step 0 complete: raw read symlinks ready ==="
