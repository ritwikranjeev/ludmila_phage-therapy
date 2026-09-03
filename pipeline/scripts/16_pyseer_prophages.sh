#!/usr/bin/env bash
set -euo pipefail

# Run pyseer on prophage features.
# Submit: sbatch slurm/16_pyseer_prophages.slurm

# Environment: slurm/*.slurm sources config/hpc_env.sh

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PHENO="${PROJECT_ROOT}/config/phenotypes.tsv"
OUT_DIR="${PROJECT_ROOT}/outputs/gwas/pyseer_prophage"
FEATURES_WIDE="${PROJECT_ROOT}/outputs/prophages/prophage_feature_matrix.tsv"
RTAB="${OUT_DIR}/prophage_features_pyseer.Rtab"
# shellcheck source=scripts/ensure_population_dist.sh
source "${PROJECT_ROOT}/scripts/ensure_population_dist.sh"

mkdir -p "${OUT_DIR}"
python3 "${PROJECT_ROOT}/scripts/build_prophage_feature_matrix.py"
[[ -s "${FEATURES_WIDE}" ]] || {
  echo "Failed to build ${FEATURES_WIDE}" >&2
  exit 1
}
python3 "${PROJECT_ROOT}/scripts/wide_presence_to_rtab.py" "${FEATURES_WIDE}" > "${RTAB}"

ensure_population_dist

for trait in binary infection_class strong_infect; do
  pyseer \
    --phenotypes "${PHENO}" \
    --phenotype-column "${trait}" \
    --pres "${RTAB}" \
    --distances "${PYSEER_STRUCTURE_DIST}" \
    > "${OUT_DIR}/pyseer_prophage_${trait}.txt"
done
