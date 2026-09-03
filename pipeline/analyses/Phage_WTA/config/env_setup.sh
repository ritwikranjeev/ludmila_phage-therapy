#!/usr/bin/env bash
# Sourced by Slurm scripts and interactive runs. Defines project paths and activates conda.
# Adjust module name or fallback paths if your HPC3 setup differs.
#
# For path-only bootstrap (no conda), e.g. bin/sanity_check_paths.sh from CI:
#   export PHAGE_WTA_SKIP_CONDA=1
# before sourcing this file.
#
# On HPC3, PHAGE_WTA_SKIP_CONDA=1 also loads an RCIC anaconda module (for python3 on
# compute nodes) and prepends HMMER on PATH — use for jobs that only need stdlib +
# external tools when phage_wta_env lives on a filesystem not mounted on workers.

PHAGE_WTA_PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PHAGE_WTA_PROJECT_ROOT
export PHAGE_WTA_DATA_DIR="${PHAGE_WTA_PROJECT_ROOT}/data"

_phage_wta_load_rcic_anaconda_module() {
  if ! command -v module &>/dev/null; then
    return 1
  fi
  if [[ -n "${PHAGE_WTA_CONDA_MODULE:-}" ]]; then
    module load "${PHAGE_WTA_CONDA_MODULE}" 2>/dev/null
    return
  fi
  for _phage_wta_mod in anaconda/2024.06 anaconda/2025.12 anaconda/2022.05 anaconda/2021.11; do
    if module load "${_phage_wta_mod}" 2>/dev/null; then
      unset _phage_wta_mod
      return 0
    fi
  done
  unset _phage_wta_mod
  return 1
}

_phage_wta_prepend_hmmer_bin() {
  if ! command -v hmmsearch &>/dev/null; then
    if [[ -n "${PHAGE_WTA_HMMER_BIN_DIR:-}" && -x "${PHAGE_WTA_HMMER_BIN_DIR}/hmmsearch" ]]; then
      export PATH="${PHAGE_WTA_HMMER_BIN_DIR}:${PATH}"
    else
      for _phage_wta_d in \
        "${HOME}/.local/share/mamba/bin" \
        "${HOME}/miniforge3/bin" \
        "${HOME}/mambaforge/bin" \
        "${HOME}/micromamba/bin"
      do
        if [[ -x "${_phage_wta_d}/hmmsearch" ]]; then
          export PATH="${_phage_wta_d}:${PATH}"
          break
        fi
      done
    fi
    unset _phage_wta_d
  fi
  if ! command -v hmmsearch &>/dev/null && command -v module &>/dev/null; then
    module load hmmer/3.3 2>/dev/null || true
  fi
}

if [[ "${PHAGE_WTA_SKIP_CONDA:-0}" == "1" ]]; then
  if ! _phage_wta_load_rcic_anaconda_module; then
    echo "WARNING: No RCIC anaconda module loaded (tried anaconda/2024.06, 2025.12, 2022.05, 2021.11). Set PHAGE_WTA_CONDA_MODULE." >&2
  fi
  _phage_wta_rcic_py_dir=
  if command -v python3 &>/dev/null; then
    _phage_wta_rcic_py_dir="$(dirname "$(command -v python3)")"
  fi
  _phage_wta_prepend_hmmer_bin
  if [[ -n "${_phage_wta_rcic_py_dir:-}" ]]; then
    export PATH="${_phage_wta_rcic_py_dir}:${PATH}"
  fi
  unset _phage_wta_rcic_py_dir
  if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    exit 0
  fi
  return 0
fi

# --- Conda: UCI HPC3 / RCIC Environment Modules ---
# RCIC LANGUAGES includes anaconda/2020.07, 2021.11, 2022.05, 2024.06, 2025.12 (see: module avail anaconda).
# Do not load other Python/tool modules in the same job as anaconda (RCIC guidance).
# Override with: export PHAGE_WTA_CONDA_MODULE="anaconda/2022.05"
if ! _phage_wta_load_rcic_anaconda_module; then
  echo "WARNING: No RCIC anaconda module loaded (tried anaconda/2024.06, 2025.12, 2022.05, 2021.11). Set PHAGE_WTA_CONDA_MODULE or use miniforge under \$HOME." >&2
fi

if command -v conda &>/dev/null; then
  # shellcheck source=/dev/null
  source "$(conda info --base)/etc/profile.d/conda.sh"
else
  for _phage_wta_conda_sh in \
    "${HOME}/miniforge3/etc/profile.d/conda.sh" \
    "${HOME}/miniconda3/etc/profile.d/conda.sh"
  do
    if [[ -f "${_phage_wta_conda_sh}" ]]; then
      # shellcheck source=/dev/null
      source "${_phage_wta_conda_sh}"
      break
    fi
  done
  unset _phage_wta_conda_sh
fi

if ! conda activate phage_wta_env; then
  {
    echo "ERROR: conda activate phage_wta_env failed."
    echo "  - Create the env: conda create -n phage_wta_env python=3.11  (on an interactive node: srun --pty bash, per RCIC)"
    echo "  - List envs: conda env list"
    echo "  - Path-only scripts: export PHAGE_WTA_SKIP_CONDA=1 before sourcing this file."
  } >&2
  return 1 2>/dev/null || exit 1
fi

_phage_wta_prepend_hmmer_bin

# Prepending a user mamba base for hmmsearch can put an old python3 first; keep env interpreter primary.
if [[ -n "${CONDA_PREFIX:-}" && -d "${CONDA_PREFIX}/bin" ]]; then
  export PATH="${CONDA_PREFIX}/bin:${PATH}"
fi
