#!/bin/bash
#SBATCH --job-name=stk2_panaroo_nr
#SBATCH --output=logs/stk2_linked_blast_%j.out
#SBATCH --error=logs/stk2_linked_blast_%j.err
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --partition=free

# Uses NCBI blastp -remote -db nr (needs outbound HTTPS from compute node).
# Submit from project root:
#   sbatch slurm/stk2_linked_panaroo_blast_remote.slurm

set -euo pipefail
PROJECT_ROOT="${PROJECT_ROOT:-/pub/ritwik1/manuscripts/OCPT-Ludmila_Characterization}"
cd "$PROJECT_ROOT"
mkdir -p logs results/reports

module purge 2>/dev/null || true
module load ncbi-blast/2.13.0 2>/dev/null || module load ncbi-blast 2>/dev/null || true
module load miniforge3 2>/dev/null || module load anaconda3 2>/dev/null || true

echo "=== Job $SLURM_JOB_ID $(hostname) ==="
date
command -v blastp
command -v python3
python3 -c "import pandas; print('pandas ok')"

GENE_DATA="$PROJECT_ROOT/outputs/pangenome/gene_data.csv"
if [[ ! -f "$GENE_DATA" ]]; then
  echo "ERROR: missing $GENE_DATA"
  exit 1
fi
if [[ ! -f "$PROJECT_ROOT/results/stk2_panaroo_gene_linkage.tsv" ]]; then
  echo "ERROR: missing results/stk2_panaroo_gene_linkage.tsv (rsync from laptop or run stk2_panaroo_gene_linkage.py)"
  exit 1
fi

# Per-sequence remote BLAST with pauses avoids NCBI SIGXCPU on multi-query jobs.
python3 "$PROJECT_ROOT/scripts/annotate_stk2_linked_panaroo.py" \
  --project-root "$PROJECT_ROOT" \
  --gene-data-csv "$GENE_DATA" \
  --blast-remote-nr \
  --remote-blast-sleep 10 \
  --remote-blast-retries 5 \
  --top-n 15

echo "=== Done ==="
date
ls -la "$PROJECT_ROOT/results/reports/stk2_linked_panaroo_annotations.tsv" || true
