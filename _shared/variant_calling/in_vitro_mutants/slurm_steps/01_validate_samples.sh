#!/bin/bash
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH -t 00:30:00
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -J fig5_validate
#SBATCH -o logs/01_validate_%j.out
#SBATCH -e logs/01_validate_%j.err
#
# Step 1: Check sample sheet and that all read paths exist under project_dir.
# Run from pipeline root, or submit with: sbatch -D /path/to/Figure5_Phage_resistant_Mutant_Analysis slurm_steps/01_validate_samples.sh

set -euo pipefail
ANALYSIS_DIR="/pub/ritwik1/manuscripts/OCPT-Ludmila_Characterization/analyses/Figure5_Phage_resistant_Mutant_Analysis"
cd "$ANALYSIS_DIR"
mkdir -p logs

module load python/3.10.2 2>/dev/null || module load python3 2>/dev/null || true

echo "=== Figure5 step 1: validate_samples.py ==="
echo "CWD: $(pwd)"
python3 scripts/validate_samples.py metadata/samples.tsv --project-dir "$ANALYSIS_DIR"
echo "=== OK: sample sheet and paths valid ==="
