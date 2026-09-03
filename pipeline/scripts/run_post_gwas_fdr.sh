#!/usr/bin/env bash
# Add BH-FDR and sort all pyseer outputs. Run from project root:
#   bash scripts/run_post_gwas_fdr.sh
set -euo pipefail
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="${PROJECT_ROOT}/scripts/18_post_gwas_stats.py"
OUT="${PROJECT_ROOT}/outputs/gwas/post_gwas_fdr"
TOP="${POST_GWAS_TOP:-50}"

mkdir -p "${OUT}"

run_one() {
  local rel="$1"
  local tag="$2"
  local in="${PROJECT_ROOT}/${rel}"
  [[ -s "${in}" ]] || {
    echo "skip (missing): ${rel}" >&2
    return 0
  }
  local base
  base="$(basename "${rel}" .txt)"
  python3 "${PY}" --infile "${in}" --outfile "${OUT}/${tag}__${base}_fdr.tsv" --top "${TOP}"
}

# gene GWAS
for t in binary infection_class strong_infect; do
  run_one "outputs/gwas/pyseer_gene/pyseer_${t}.txt" "gene_${t}"
done
# defense
for t in binary infection_class strong_infect; do
  run_one "outputs/gwas/pyseer_defense/pyseer_defense_${t}.txt" "defense_${t}"
done
# prophage
for t in binary infection_class strong_infect; do
  run_one "outputs/gwas/pyseer_prophage/pyseer_prophage_${t}.txt" "prophage_${t}"
done
# host receptors
for t in binary infection_class strong_infect; do
  run_one "outputs/gwas/pyseer_host_receptors/pyseer_host_receptors_${t}.txt" "hostrec_${t}"
done

python3 "${PROJECT_ROOT}/scripts/merge_gwas_priority_overview.py" \
  --indir "${OUT}" \
  --outfile "${PROJECT_ROOT}/outputs/gwas/gwas_priority_overview.tsv" \
  --per-file 15

echo "Done. FDR tables under ${OUT}/ (top ${TOP} per analysis; set POST_GWAS_TOP to change)."
echo "Merged preview: outputs/gwas/gwas_priority_overview.tsv"
