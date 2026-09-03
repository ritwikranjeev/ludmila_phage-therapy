# Copy or symlink to hpc_env.sh only if you want a local override.
# Committed defaults: config/hpc_env.sh (edit CONDA_ENV_* there or via export before sbatch).

export MINIFORGE_ROOT="/pub/ritwik1/housekeeping/miniforge3"
conda_activate() {
  local env_name="$1"
  local bindir="${MINIFORGE_ROOT}/envs/${env_name}/bin"
  export PATH="${bindir}:${PATH}"
  export CONDA_PREFIX="${MINIFORGE_ROOT}/envs/${env_name}"
}
