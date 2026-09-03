#!/usr/bin/env bash
set -euo pipefail

# Run pyseer on defense island features.
# Submit: sbatch slurm/15_pyseer_defense_islands.slurm

# Environment: slurm/*.slurm sources config/hpc_env.sh

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PHENO="${PROJECT_ROOT}/config/phenotypes.tsv"
OUT_DIR="${PROJECT_ROOT}/outputs/gwas/pyseer_defense"
FEATURES_WIDE="${PROJECT_ROOT}/outputs/defense/defense_family_presence.tsv"
# pyseer --pres expects roary/panaroo Rtab: rows = features, columns = samples.
RTAB="${OUT_DIR}/defense_family_presence_pyseer.Rtab"
# shellcheck source=scripts/ensure_population_dist.sh
source "${PROJECT_ROOT}/scripts/ensure_population_dist.sh"

mkdir -p "${OUT_DIR}"
[[ -s "${FEATURES_WIDE}" ]] || {
  echo "Missing ${FEATURES_WIDE} (run step 10 defense islands first)" >&2
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
    > "${OUT_DIR}/pyseer_defense_${trait}.txt"
done
