# Figure 2a Ludmila phage assembly pipeline – shared config
# Source this from scripts: source "$(dirname "$0")/../config/config.sh"
# Or from pipeline root: source config/config.sh

CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PROJECT_ROOT="$(cd "$CONFIG_DIR/.." && pwd)"

# Phage name (used for read naming, assembly output, etc.)
export PHAGE_NAME="${PHAGE_NAME:-Ludmila}"

# Data and results (under project root)
export DATA_DIR="$PROJECT_ROOT/data"
export RAW_READS_DIR="$DATA_DIR/raw_fastqs"
export RESULTS_DIR="$PROJECT_ROOT/results"

# Cleaned reads and QC
export READS_CLEAN="$RESULTS_DIR/reads/clean"
export QC_DIR="$RESULTS_DIR/qc"

# Assembly
export ASSEMBLY_DIR="$RESULTS_DIR/assembly"
export ASSEMBLY_FASTA="$ASSEMBLY_DIR/${PHAGE_NAME}.fasta"

# CheckV
export CHECKV_DIR="$RESULTS_DIR/checkv"
export CHECKV_DB="${CHECKV_DB:-/dfs6b/pub/ritwik1/housekeeping/db/checkv/checkv-db-v1.5}"

# Annotation (Pharokka → Phold)
export ANNOT_DIR="$RESULTS_DIR/annotation"
export PHAROKKA_DIR="$ANNOT_DIR/pharokka"
export PHAROKKA_GBK="$PHAROKKA_DIR/pharokka.gbk"   # Pharokka output GenBank
export PHOLD_DIR="$ANNOT_DIR/phold"
export PHOLD_COMPARE_DIR="$ANNOT_DIR/phold/compare"
export PHOLD_PLOT_DIR="$ANNOT_DIR/phold/plot"

# Optional: PhiX and host reference for bbduk (leave empty to skip)
export PHIX_FASTA="${PHIX_FASTA:-}"
export HOST_FASTA="${HOST_FASTA:-}"

# Pharokka database (required for pharokka_annotation.sh).
# Pharokka 1.9+ expects v1.8.0 DB (VERSION_1_8_0, INPHARED mash files). If the default path
# has an older DB, run: install_databases.py -o /path/to/dir  (from pharokkaENV), then set
# PHAROKKA_DB to that path, or set PHAROKKA_DB to an existing v1.8 DB directory.
export PHAROKKA_DB="${PHAROKKA_DB:-/dfs6b/pub/ritwik1/housekeeping/db/pharokka_v1.8.0_databases}"

# Phold plot title (optional)
export PHOLD_PLOT_TITLE="${PHOLD_PLOT_TITLE:-S. aureus Phage Ludmila24}"

# Optional: path to quast.py for assembly QC (leave empty to skip QUAST)
export QUAST_PYTHON="${QUAST_PYTHON:-}"

# Optional: conda init script for seqkit in assembly.sh (e.g. /path/to/miniforge3/etc/profile.d/conda.sh)
export CONDA_INIT="${CONDA_INIT:-}"

# Miniforge on HPC3. This install was moved under /dfs6b/ but many bin/* shebangs still point at
# /dfs6/... (bad interpreter). Do not rely on conda activate or #!/.../bin/phold; call each
# interpreter explicitly (see PHOLD_PYTHON + PHOLD_BIN, etc.).
export MINIFORGE_ROOT="${MINIFORGE_ROOT:-/dfs6b/pub/ritwik1/housekeeping/miniforge3}"

export PHOLD_PYTHON="${PHOLD_PYTHON:-${MINIFORGE_ROOT}/envs/pholdENV/bin/python3.11}"
export PHOLD_BIN="${PHOLD_BIN:-${MINIFORGE_ROOT}/envs/pholdENV/bin/phold}"
export PHAROKKA_PYTHON="${PHAROKKA_PYTHON:-${MINIFORGE_ROOT}/envs/pharokkaENV/bin/python3.13}"
export PHAROKKA_BIN="${PHAROKKA_BIN:-${MINIFORGE_ROOT}/envs/pharokkaENV/bin/pharokka.py}"
export CHECKV_PYTHON="${CHECKV_PYTHON:-${MINIFORGE_ROOT}/envs/checkV/bin/python3.12}"
export CHECKV_BIN="${CHECKV_BIN:-${MINIFORGE_ROOT}/envs/checkV/bin/checkv}"

# SLURM (edit for your account)
export SLURM_ACCOUNT="${SLURM_ACCOUNT:-katrine_lab}"
export SLURM_PARTITION="${SLURM_PARTITION:-standard}"
export SLURM_MAIL_USER="${SLURM_MAIL_USER:-ritwik1@uci.edu}"
