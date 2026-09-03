#!/bin/bash
#SBATCH --job-name=stk2_panaroo_dia
#SBATCH --output=logs/stk2_linked_diamond_%j.out
#SBATCH --error=logs/stk2_linked_diamond_%j.err
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
## HPC3 (RCIC): pick a partition you can use, e.g. free, debug, or your lab QOS
#SBATCH --partition=free

# Stk2-linked Panaroo genes: DIAMOND blastp + build annotation TSVs.
#
# Submit from the project root (or set PROJECT_ROOT):
#   export DIAMOND_DB=/path/to/uniref50.dmnd    # required
#   export GENE_DATA=/path/to/gene_data.csv     # optional if file is under outputs/pangenome/
#   sbatch slurm/stk2_linked_panaroo_diamond.slurm
#
# Or one-liner:
#   sbatch --export=ALL,DIAMOND_DB=$HOME/db/uniref50.dmnd slurm/stk2_linked_panaroo_diamond.slurm

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"
cd "$PROJECT_ROOT"

mkdir -p logs results/reports

echo "=== Job $SLURM_JOB_ID node $(hostname) ==="
echo "PROJECT_ROOT=$PROJECT_ROOT"
date

if [[ -z "${DIAMOND_DB:-}" ]]; then
  echo "ERROR: Export DIAMOND_DB=/path/to/database.dmnd before sbatch."
  exit 1
fi
if [[ ! -f "$DIAMOND_DB" ]]; then
  echo "ERROR: DIAMOND_DB file not found: $DIAMOND_DB"
  exit 1
fi

GENE_DATA="${GENE_DATA:-$PROJECT_ROOT/outputs/pangenome/gene_data.csv}"
COMBINED_FA="${COMBINED_FA:-$PROJECT_ROOT/outputs/pangenome/combined_protein_CDS.fasta}"

EXTRA_ARGS=()
if [[ -f "$GENE_DATA" ]]; then
  EXTRA_ARGS+=(--gene-data-csv "$GENE_DATA")
  echo "Using gene_data: $GENE_DATA"
elif [[ -f "$COMBINED_FA" ]]; then
  EXTRA_ARGS+=(--combined-protein-fasta "$COMBINED_FA")
  echo "Using combined FASTA: $COMBINED_FA"
else
  echo "ERROR: Need Panaroo gene_data.csv or combined_protein_CDS.fasta under outputs/pangenome/, or set GENE_DATA / COMBINED_FA."
  exit 1
fi

# --- Modules: adjust for your HPC3 stack ---
module purge 2>/dev/null || true
module load diamond 2>/dev/null || module load DIAMOND 2>/dev/null || true
# Python with pandas + statsmodels if needed (annotate script only needs pandas)
module load miniforge3 2>/dev/null || module load anaconda3 2>/dev/null || true
if command -v conda >/dev/null 2>&1; then
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null && conda activate base 2>/dev/null || true
fi

THR="${SLURM_CPUS_PER_TASK:-8}"
export PATH="${PATH:-}"

if ! command -v diamond >/dev/null 2>&1; then
  echo "ERROR: diamond not on PATH after modules. Load the correct module or install DIAMOND in your env."
  exit 1
fi
if ! command -v python3 >/dev/null 2>&1; then
  echo "ERROR: python3 not on PATH."
  exit 1
fi

diamond version
python3 -c "import pandas; print('pandas OK')"

python3 "$PROJECT_ROOT/scripts/annotate_stk2_linked_panaroo.py" \
  --project-root "$PROJECT_ROOT" \
  "${EXTRA_ARGS[@]}" \
  --diamond-db "$DIAMOND_DB" \
  --threads "$THR" \
  --top-n "${TOP_N:-15}"

echo "=== Done ==="
date
ls -la "$PROJECT_ROOT/results/reports/stk2_linked_panaroo_annotations.tsv" 2>/dev/null || true
