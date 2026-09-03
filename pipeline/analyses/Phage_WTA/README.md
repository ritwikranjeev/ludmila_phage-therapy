# Phage–WTA analysis (Ludmila)

This directory holds a minimal scaffold for phage–wall teichoic acid (WTA) analysis work.

**Step-by-step HPC3 runbook:** see **[PIPELINE_STEPS.md](PIPELINE_STEPS.md)** for the full dry-lab order and concrete commands. Run **`bash bin/sanity_check_paths.sh`** once to ensure expected **`logs/`**, **`data/`**, and **`analysis/`** directories exist. On the cluster, **`bash bin/print_pipeline_status.sh`** prints **FOUND**/**MISSING** for key outputs (it does not submit jobs). Tiny **dummy** phage/host FASTAs (`data/phage/Ludmila_test.fa`, `data/hosts/HostA_test.fa`) and matching **`config/*_genomes.tsv`** rows support wiring checks only—replace with real assemblies before science runs.

**Before first real run (checklist):**

1. Confirm **`config/env_setup.sh`** loads the right **conda** / **module** stack and that **`phage_wta_env`** exists (`conda activate phage_wta_env` works on the login node).
2. Replace **MTESTSEQ** placeholders in **`config/wta_reference_genes.fasta`** with curated references (strict **`python bin/preflight_check.py`** fails until you do).
3. Set **GPU/CUDA module** lines in **`slurm/run_colabfold_ludmila_rbp.sh`** per RCIC docs. ColabFold jobs target **`free-gpu`** by default (no GPU lab account; preemptible).
4. Verify **`colabfold_batch`** and **`prokka`** (and **HMMER** binaries) are on **PATH** inside **`phage_wta_env`**.
5. Run **`export DRY_RUN=1`** and a **single-array** **Prokka** job (and optionally **`sbatch slurm/preflight_check.sh`**); inspect **`logs/`** for echoed commands, then **`unset DRY_RUN`**. **`slurm/template_cpu_job.sh`** / **`template_gpu_job.sh`** also honor **`DRY_RUN=1`** for copy-paste tests.

**Quick checks:** **`python3 bin/preflight_check.py`** (add **`--allow-wta-placeholder`** only while references are still dummy), **`bash bin/run_smoke_test.sh`** (local wiring, no **Slurm**), **`python3 bin/validate_troubleshoot_complete.py`** (full local gate; **`--print-hpc3-commands`** for SSH/sbatch text). Override the printed SSH path with **`PHAGE_WTA_HPC3_USER`**, **`PHAGE_WTA_HPC3_HOST`**, and **`PHAGE_WTA_HPC3_PROJECT_DIR`** (default remote dir: **`/pub/ritwik1/manuscripts/OCPT-Ludmila_PhageTherapy/WTA-RBP_Analysis`**). Audit notes: **`analysis/repo_audit.md`**.

**Pipeline stages (0–6):** setup and preflight → genome lists → Prokka + optional RBP tool → WTA HMM build + host scan + feature matrix → Ludmila RBP extraction + ColabFold → infectivity association → inspect outputs. Step-by-step commands live in **`PIPELINE_STEPS.md`**.

**`raw_wta_hits.tsv`** columns (from **`scan_hosts_for_wta_genes.py`**): **`host_id`**, **`gene`**, **`best_evalue`**, **`bitscore`**, **`hit_length`**. Re-run scans if you still have older **`by_host/*_wta_hits.tsv`** files with the previous column order.

**Environment hints (HPC3):** RCIC **`LANGUAGES`** includes **`anaconda/2024.06`**, **`anaconda/2025.12`**, **`anaconda/2022.05`**, **`anaconda/2021.11`**, **`anaconda/2020.07`**, plus **`bioconda`**, **`mamba`**, **`miniconda3`**. Override with **`export PHAGE_WTA_CONDA_MODULE=anaconda/2025.12`** (or another) before jobs if the default order in **`env_setup.sh`** is wrong. **HMMER** and **BLAST+** exist as **`BIOTOOLS`** modules (**`hmmer/3.3`**, **`ncbi-blast/2.13.0`**) if you prefer cluster installs over conda—see **`PIPELINE_STEPS.md`** for PATH caveats. **Prokka** is not a standalone RCIC module in the usual listing; use **Bioconda** in **`phage_wta_env`**. For **`associate_wta_with_infectivity.py`** on very small tables, set **`export PHAGE_WTA_ASSOC_MIN_HOSTS=2`** or **`python bin/associate_wta_with_infectivity.py --min-hosts 2`**.

**Conda install sketch** (after **`module load anaconda/2024.06`** (or another RCIC **`anaconda/...`**) and **`conda create -n phage_wta_env python=3.11`**, on an interactive node per RCIC): activate the env, then install **Prokka** / **HMMER** / **blast** from **Bioconda** as needed, **`pip install pandas scikit-learn numpy`**, and **ColabFold** exactly as required by the ColabFold release you choose (verify GPU/JAX/CUDA compatibility on HPC3).

Jobs are intended to run on **UCI HPC3** with **Slurm**, using Slurm account **`katrine_lab`**.

## First steps after cloning

1. Create a conda environment named **`phage_wta_env`** (tools and versions are up to you; this repo does not create the env).
2. Open `config/env_setup.sh` and adjust conda loading if needed (module name, or rely on a personal Miniforge/Miniconda install under `$HOME`).
3. Place or sync genomes under `data/` as needed.

## Step 1: Annotation and RBP detection

- Put phage assemblies in `data/phage/` as `.fasta` or `.fa`, and host assemblies in `data/hosts/` with the same extensions.
- From the project root, generate Slurm input lists:
  ```bash
  python bin/prepare_genome_lists.py
  ```
  This writes `config/phage_genomes.tsv` and `config/host_genomes.tsv` (columns `id`, `path`). Missing `data/phage/` or `data/hosts/` only produces a warning.
- Submit **array** jobs (set the upper bound from the number of **data** lines in each TSV, i.e. total lines minus one for the header). Examples — **replace `N` and `M` with your counts** (e.g. `N=$(($(wc -l < config/phage_genomes.tsv) - 1))` then use `1-$N`):
  ```bash
  # Phage Prokka (one task per phage genome)
  sbatch --array=1-N slurm/annotate_phage_prokka.sh

  # Host Prokka (one task per host genome)
  sbatch --array=1-M slurm/annotate_hosts_prokka.sh

  # Phage RBP detection (requires Prokka outputs under data/phage_prokka/<id>/; run after phage annotation)
  sbatch --array=1-N slurm/detect_phage_rbp.sh
  ```
- Install **Prokka** and **PhageRBPdetection** in `phage_wta_env` (or module PATH). Edit `slurm/detect_phage_rbp.sh` TODO block once you confirm the real PhageRBPdetection command line.

## Step 2: WTA feature inference

- Replace dummy sequences in `config/wta_reference_genes.fasta` with curated WTA pathway references, then build the HMM library once on the cluster:
  ```bash
  sbatch slurm/build_wta_hmmdb.sh
  ```
  This runs `python bin/build_wta_hmmdb.py`, which writes per-gene HMMs and a combined database under `data/wta_hmms/` plus `hmmpress` index files (`.h3f`, `.h3i`, `.h3m`, `.h3p`).
- Complete **Step 1** host annotation first so each host has `data/host_prokka/<id>/<id>.faa` (or adjust paths in `bin/scan_hosts_for_wta_genes.py`).
- Scan hosts in parallel (**replace `M` with your host count**, same rule as Step 1: number of data lines in `config/host_genomes.tsv`):
  ```bash
  sbatch --array=1-M slurm/scan_hosts_for_wta_genes.sh
  ```
  Each task writes `data/wta_features/by_host/<id>_wta_hits.tsv` and `data/wta_features/by_host/<id>.domtblout` (combined-HMM mode).
- Merge array outputs, then build the matrix:
  ```bash
  python bin/scan_hosts_for_wta_genes.py --collect-raw
  python bin/build_wta_feature_matrix.py
  ```
  The matrix is `data/wta_features/wta_feature_matrix.tsv`. Refine gene→backbone/glycan rules in `build_wta_feature_matrix.py` where marked **TODO**.

## Step 3: RBP structure and WTA–infectivity association

- After Ludmila (or your phage of interest) is annotated under `data/phage_prokka/<name>/`, build candidate tail/RBP FASTA for ColabFold:
  ```bash
  python bin/extract_ludmila_rbp_candidates.py \
    --prokka-dir data/phage_prokka/Ludmila \
    --out data/phage_rbp/Ludmila_RBP_candidates.fasta
  ```
  Adjust `--prokka-dir` / `--out` and, if needed, protein filename (e.g. `Ludmila.faa`) by keeping the Prokka prefix aligned with the directory name.
- Submit GPU structure prediction (edit CUDA modules and ColabFold flags in the script as needed):
  ```bash
  sbatch slurm/run_colabfold_ludmila_rbp.sh
  ```
  Models are written under `data/alphafold/Ludmila_RBPs/`.
- Prepare infectivity labels (creates a template if missing, or validates an existing file):
  ```bash
  python bin/prepare_infectivity_table.py
  ```
  Replace dummy rows in `data/phenotypes/ludmila_infectivity.tsv` with real Ludmila infectivity data; `host_id` must match `wta_feature_matrix.tsv`.
- Run association analysis (CPU job or locally):
  ```bash
  sbatch slurm/associate_wta_with_infectivity.sh
  ```
  Outputs `analysis/wta_infectivity_logreg_coeffs.tsv` and prints accuracy, ROC AUC (when defined), and coefficient ranking.

## Layout

| Path        | Purpose                                      |
|------------|-----------------------------------------------|
| `bin/`     | Small helper scripts (e.g. `prepare_genome_lists.py`) |
| `slurm/`   | Slurm batch templates and array drivers       |
| `config/`  | Paths, `env_setup.sh`, genome TSV lists       |
| `analysis/`| Notebooks or downstream analysis code         |
| `data/`    | `phage/`, `hosts/`; annotations `phage_prokka/`, `host_prokka/`, `phage_rbp/`; WTA `wta_hmms/`, `wta_features/`; phenotypes `phenotypes/`; structures `alphafold/` |

**Slurm billing (CPU):** every CPU driver under **`slurm/`** uses **`#SBATCH --account=katrine_lab`** and **`--partition=standard`** (Prokka, WTA HMM build/scan, association, preflight). **GPU** ColabFold uses **`free-gpu`**; it also sets **`--account=katrine_lab`** for now (not charged like **`gpu`**). Override at submit time with **`sbatch --account=...`** only if your PI asks.

Submit jobs from this project root so Slurm sets **`SLURM_SUBMIT_DIR`** correctly (batch scripts **`cd`** there before sourcing **`config/env_setup.sh`**). Example:

```bash
cd /pub/ritwik1/manuscripts/OCPT-Ludmila_PhageTherapy/WTA-RBP_Analysis
sbatch slurm/template_cpu_job.sh
```
