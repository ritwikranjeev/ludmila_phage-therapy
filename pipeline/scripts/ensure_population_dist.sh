#!/usr/bin/env bash
# Source from pyseer bash scripts after PROJECT_ROOT and PHENO are set:
#   # shellcheck source=scripts/ensure_population_dist.sh
#   source "${PROJECT_ROOT}/scripts/ensure_population_dist.sh"
#   ensure_population_dist

# Writes outputs/gwas/population_structure_dist.tsv from IQ-TREE Newick (patristic distances).
# Requires: pip install dendropy

ensure_population_dist() {
  local tree="${PROJECT_ROOT}/outputs/trees/core.treefile"
  local struct="${PROJECT_ROOT}/outputs/gwas/population_structure_dist.tsv"
  [[ -n "${PHENO:-}" ]] || {
    echo "ensure_population_dist: set PHENO to config/phenotypes.tsv" >&2
    return 1
  }
  [[ -s "${tree}" ]] || {
    echo "Missing IQ-TREE tree: ${tree}" >&2
    return 1
  }
  mkdir -p "$(dirname "${struct}")"
  if [[ ! -s "${struct}" ]] || [[ "${tree}" -nt "${struct}" ]]; then
    python3 "${PROJECT_ROOT}/scripts/phylogeny_distance.py" \
      "${tree}" --phenotypes-tsv "${PHENO}" > "${struct}"
  fi
  export PYSEER_STRUCTURE_DIST="${struct}"
}
