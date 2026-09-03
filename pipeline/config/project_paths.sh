# =============================================================================
# project_paths.sh -- Shared project root and paths for Ludmila manuscript.
# Source this from scripts so PROJECT_DIR points at LudmilaCharacterization/.
# =============================================================================

# PROJECT_DIR: manuscript project root (LudmilaCharacterization/)
# If not set, infer from this script's location (config/project_paths.sh -> parent = project root).
if [[ -z "${PROJECT_DIR:-}" ]]; then
    _SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-.}")" && pwd)"
    PROJECT_DIR="$(cd "$_SCRIPT_DIR/.." && pwd)"
    unset _SCRIPT_DIR
fi
export PROJECT_DIR

# Central data and analysis roots (all under PROJECT_DIR)
DATA_DIR="${DATA_DIR:-$PROJECT_DIR/data}"
RESULTS_DIR="${RESULTS_DIR:-$PROJECT_DIR/analyses}"

# Optional HPC/cluster overrides (set in pipeline_config.sh or env)
CONDA_INIT_SCRIPT="${CONDA_INIT_SCRIPT:-}"
CONDA_ENV="${CONDA_ENV:-}"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-}"

# For host-range (and other analyses), the driver sets:
#   DATA_DIR=$PROJECT_DIR/data/phage_hostrange
#   RESULTS_DIR=$PROJECT_DIR/analyses/phage_hostrange/results
export DATA_DIR
export RESULTS_DIR
