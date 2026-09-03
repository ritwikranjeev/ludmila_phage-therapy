# Source from bash: set +e; source scripts/activate_conda_snakemake.sh; set -e
# Snakemake --use-conda shells out to conda; ensure bootstrap install is on PATH + CONDA_EXE.
find_conda_sh_fig5() {
  local c
  for c in \
    "$HOME/miniforge3/etc/profile.d/conda.sh" \
    "$HOME/miniconda3_snakemake/etc/profile.d/conda.sh" \
    "$HOME/miniconda3/etc/profile.d/conda.sh"
  do
    [[ -f "$c" ]] || continue
    echo "$c"
    return 0
  done
  return 1
}

if ! CONDA_SH_FIG5="$(find_conda_sh_fig5)"; then
  echo "ERROR: No usable conda profile.d/conda.sh under HOME." >&2
  return 1
fi

CONDA_ROOT="$(dirname "$(dirname "$CONDA_SH_FIG5")")"
export CONDA_ROOT
export CONDA_EXE="${CONDA_EXE:-$CONDA_ROOT/bin/conda}"
export PATH="$CONDA_ROOT/bin:$PATH"

SNAKEMAKE_ENV="${SNAKEMAKE_CONDA_ENV:-/dfs6b/pub/ritwik1/housekeeping/miniforge3/envs/snakemake}"

# shellcheck source=/dev/null
source "$CONDA_SH_FIG5"
conda activate "$SNAKEMAKE_ENV"
