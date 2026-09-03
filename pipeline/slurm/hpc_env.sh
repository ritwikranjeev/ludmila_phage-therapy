#!/usr/bin/env bash
# HPC3 defaults for ritwik1 — source from SLURM wrappers or interactively:
#   source /pub/ritwik1/manuscripts/OCPT-Ludmila_Characterization/config/hpc_env.sh
#   hpc_setup_step_05_prokka
#
# Override conda env names if an env is broken, e.g.:
#   export CONDA_ENV_DEFENSE=defensefinderENV

export MINIFORGE_ROOT="${MINIFORGE_ROOT:-/pub/ritwik1/housekeeping/miniforge3}"

conda_activate() {
  local env_name="$1"
  local root="${MINIFORGE_ROOT}/envs/${env_name}"
  local bindir="${root}/bin"
  if [[ ! -d "${bindir}" ]]; then
    echo "conda_activate: missing env directory: ${root}" >&2
    return 1
  fi
  # Avoid conda/conda.sh here: many HPC installs hardcode a dead /dfs6/... prefix in those scripts.
  export PATH="${bindir}:${PATH}"
  export CONDA_PREFIX="${root}"
}

# --- Conda env names (from your envs/ list; swap if a job fails) ---
export CONDA_ENV_LONG_QC_ASM="${CONDA_ENV_LONG_QC_ASM:-nanopore_env}"
export CONDA_ENV_PROKKA="${CONDA_ENV_PROKKA:-prokka}"
export CONDA_ENV_QUAST="${CONDA_ENV_QUAST:-quast}"
export CONDA_ENV_CHECKM="${CONDA_ENV_CHECKM:-checkm}"
export CONDA_ENV_PANAROO="${CONDA_ENV_PANAROO:-panarooENV}"
export CONDA_ENV_IQTREE="${CONDA_ENV_IQTREE:-iqtreeENV}"
export CONDA_ENV_DEFENSE="${CONDA_ENV_DEFENSE:-defensefinder_conda}"
export CONDA_ENV_PROPHAGE="${CONDA_ENV_PROPHAGE:-phage_defense}"
export CONDA_ENV_PYSEER="${CONDA_ENV_PYSEER:-saureus_gwas}"

# 01 — fetch/copies only (rsync/cp from login path)
hpc_setup_step_01_fetch() { :; }

# 02 — FastQC + fastp (modules). Long reads: FastQC in-script; optional NanoPlot if on PATH.
# nanopore_env has flye but may lack NanoPlot — QC does not depend on conda here.
hpc_setup_step_02_qc() {
  module load fastqc/0.11.9
  module load fastp/0.23.2
}

# 03 — SPAdes + Unicycler (modules); Flye (conda — not on module avail)
hpc_setup_step_03_assemble() {
  module load SPAdes/3.15.4
  module load unicycler/0.5.0
  module load samtools/1.15.1
  conda_activate "${CONDA_ENV_LONG_QC_ASM}"
}

# 04
hpc_setup_step_04_quast() {
  conda_activate "${CONDA_ENV_QUAST}"
}

# 04b — CheckM (Python; keep env Perl isolated from user PERL5LIB)
hpc_setup_step_04b_checkm() {
  conda_activate "${CONDA_ENV_CHECKM}"
  unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MM_OPT PERL_MB_OPT
}

# 04c — assembly cleanup (minimap2, bwa, samtools modules)
hpc_setup_step_04c_assembly_cleanup() {
  module load gcc/11.2.0
  module load minimap2/2.28
  module load bwa/0.7.17
  module load samtools/1.15.1
}

# 05a — read-level contamination diagnostics (same tools as 04c)
hpc_setup_step_05a_diagnose_contamination() {
  module load gcc/11.2.0
  module load minimap2/2.28
  module load bwa/0.7.17
  module load samtools/1.15.1
}

# 05
hpc_setup_step_05_prokka() {
  conda_activate "${CONDA_ENV_PROKKA}"
  # Login/interactive shells often set PERL5LIB to ~/perl5; that prepends broken @INC on compute nodes.
  unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MM_OPT PERL_MB_OPT
}

# 06
hpc_setup_step_06_panaroo() {
  conda_activate "${CONDA_ENV_PANAROO}"
}

# 07
hpc_setup_step_07_iqtree() {
  conda_activate "${CONDA_ENV_IQTREE}"
}

# 08 — DefenseFinder (defense-finder CLI)
hpc_setup_step_08_defense() {
  conda_activate "${CONDA_ENV_DEFENSE}"
}

# 09–10, 12–13 — stdlib Python only (uses login python3; may be 3.6)
hpc_setup_step_python_stdlib() { :; }

# 19 — host receptor scan (needs Python ≥3.7 for __future__ annotations)
hpc_setup_step_19_host_receptors() {
  module load python/3.10.2
}

# 11 — PhiSpy
hpc_setup_step_11_prophage() {
  conda_activate "${CONDA_ENV_PROPHAGE}"
}

# 14–16 — pyseer / GWAS
hpc_setup_step_pyseer() {
  conda_activate "${CONDA_ENV_PYSEER}"
}
