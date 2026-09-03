#!/bin/bash
# Local wiring smoke test: directories, genome lists, preflight (does not sbatch jobs).
# Run from project root:  bash bin/run_smoke_test.sh

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

export PHAGE_WTA_PROJECT_ROOT="${ROOT}"
export PHAGE_WTA_SKIP_CONDA=1

echo "=== 1. sanity_check_paths.sh ==="
bash "${ROOT}/bin/sanity_check_paths.sh"

echo ""
echo "=== 2. prepare_genome_lists.py ==="
python3 "${ROOT}/bin/prepare_genome_lists.py"

echo ""
echo "=== 3. preflight_check.py (allow WTA placeholder for scaffold repo) ==="
python3 "${ROOT}/bin/preflight_check.py" --allow-wta-placeholder

echo ""
echo "=== 4. build_wta_feature_matrix.py --allow-empty-raw (no HMM scan required) ==="
python3 "${ROOT}/bin/build_wta_feature_matrix.py" --allow-empty-raw

echo ""
echo "=== Suggested next steps on HPC3 (login node, project root) ==="
echo "  export DRY_RUN=1"
echo "  sbatch --array=1-1 slurm/annotate_phage_prokka.sh"
echo "  sbatch --array=1-1 slurm/annotate_hosts_prokka.sh"
echo "  # Inspect logs/phage_annotate_*_1.out and logs/host_annotate_*_1.out"
echo ""
echo "  unset DRY_RUN"
echo "  sbatch --array=1-1 slurm/annotate_phage_prokka.sh   # real Prokka (tiny test genome)"
echo "  sbatch --array=1-1 slurm/annotate_hosts_prokka.sh"
echo ""
echo "  sbatch slurm/preflight_check.sh"
echo "  # Strict preflight (fails until MTESTSEQ removed): python3 bin/preflight_check.py"
echo ""
echo "Smoke test finished (no Slurm jobs submitted from this script)."
