#!/usr/bin/env bash
# Step 1: create conda env phage_wta_env on HPC3 with CPU-pipeline dependencies.
# Run inside an interactive compute session (recommended): large conda transactions on
# shared /pub can be slow or fail on login nodes; use e.g.
#   srun -A katrine_lab -p standard --time=4:00:00 --mem=32G --cpus-per-task=4 --pty bash
# then:
#   bash bin/hpc3_create_phage_wta_env.sh
set -eo pipefail

export TMPDIR="${TMPDIR:-/tmp}"
mkdir -p "${TMPDIR}/conda_work"
export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-${TMPDIR}/conda_pkgs_phage_wta}"
mkdir -p "${CONDA_PKGS_DIRS}"

if command -v module &>/dev/null; then
  # shellcheck source=/dev/null
  source /opt/rcic/Modules/init/bash 2>/dev/null || true
  module load anaconda/2024.06
fi
# shellcheck source=/dev/null
source "$(conda info --base)/etc/profile.d/conda.sh"

ENV_NAME="${PHAGE_WTA_CONDA_ENV_NAME:-phage_wta_env}"
ENV_PREFIX="$(conda env list | awk -v n="${ENV_NAME}" '$1==n {print $NF; exit}')"
if [[ -n "${ENV_PREFIX}" ]]; then
  echo "Removing existing env: ${ENV_NAME} at ${ENV_PREFIX}"
  conda env remove -n "${ENV_NAME}" -y
fi

echo "Creating ${ENV_NAME} (python 3.11 + numpy + pandas + scikit-learn) ..."
conda create -n "${ENV_NAME}" -c conda-forge \
  python=3.11 numpy pandas scikit-learn -y

PREFIX="$(conda env list | awk -v n="${ENV_NAME}" '$1==n {print $NF; exit}')"
[[ -n "${PREFIX}" ]] || { echo "ERROR: conda env ${ENV_NAME} not found after create"; exit 1; }
[[ -x "${PREFIX}/bin/python" ]] || { echo "ERROR: ${PREFIX}/bin/python missing (partial install?)"; exit 1; }
"${PREFIX}/bin/python" -V

echo "Installing hmmer, blast, prokka (bioconda) ..."
conda install -n "${ENV_NAME}" -c conda-forge -c bioconda hmmer blast prokka -y

conda activate "${ENV_NAME}"
which python
python -V
command -v prokka >/dev/null && prokka --version | head -1
command -v hmmsearch >/dev/null && hmmsearch -h 2>&1 | head -1
command -v blastp >/dev/null && blastp -version 2>&1 | head -1
python -c "import numpy, pandas, sklearn; print('imports OK', numpy.__version__, pandas.__version__, sklearn.__version__)"
echo "Done: conda activate ${ENV_NAME}"
