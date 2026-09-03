# Phage–WTA dry-lab pipeline (HPC3)

Run these steps **on UCI HPC3** from the **Phage_WTA project directory** (the folder that contains `bin/`, `config/`, `slurm/`, and `data/`). Adjust paths if your clone lives elsewhere.

**Stage map:** **0** setup · **1** genome lists · **2** Prokka + optional phage RBP tool · **3** WTA HMMs + host scan + matrix · **4** Ludmila RBP FASTA + ColabFold · **5** infectivity association · **6** inspect outputs.

Quick progress check (no jobs submitted):

```bash
bash bin/print_pipeline_status.sh
```

---

## 0. Setup

- **Conda environment** (create once; install tools you need inside it):

  ```bash
  conda create -n phage_wta_env python=3.11 -y
  conda activate phage_wta_env
  ```

  You must install the actual dependencies yourself, for example: **Prokka**, **HMMER** (`hmmbuild`, `hmmpress`, `hmmsearch`), **pandas**, **scikit-learn**, **ColabFold** / `colabfold_batch`, plus any **PhageRBPdetection** entry point you use. Use `conda`, `pip`, or cluster modules as appropriate.

- **Project environment**: every Slurm script sources `config/env_setup.sh`, which activates `phage_wta_env` and sets `PHAGE_WTA_PROJECT_ROOT`. If HPC3 uses a different module or conda layout, **edit `config/env_setup.sh`** (module name, `conda.sh` path, etc.).

- **RCIC modules that match this pipeline** (from `module avail`; names change over time—use `module spider` if a version is retired):
  - **Python / conda:** `LANGUAGES` — **`anaconda/2024.06`**, **`anaconda/2025.12`**, **`anaconda/2022.05`**, **`anaconda/2021.11`**, **`anaconda/2020.07`**; also **`bioconda/4.8.3`**, **`mamba/24.3.0`**, **`mamba/26.1.0`**, **`miniconda3/...`**. `config/env_setup.sh` tries several **`anaconda/...`** versions in order unless **`PHAGE_WTA_CONDA_MODULE`** is set.
  - **HMMER + BLAST+ (optional alternative to conda):** `BIOTOOLS` — **`hmmer/3.3`**, **`ncbi-blast/2.13.0`** (or **`ncbi-blast/2.10.0`**). You can `module load` these in an interactive shell or add **`module load ...`** lines to Slurm scripts **if** you are not relying on the same tools from **`phage_wta_env`** (avoid duplicate/conflicting **`hmmsearch`** / **`blastp`** on **`PATH`** without checking order).
  - **Prokka:** **not** listed as its own RCIC module in the usual tree — plan on **Bioconda** / **`conda install -c bioconda prokka`** inside **`phage_wta_env`** (or another env you activate in jobs).
  - **GPU / ColabFold:** `TOOLS` — **`cuda/11.7.1`**, **`cuda/12.2.0`**, **`cuda/13.0.1`**, etc. Uncomment and pin the pair that matches your ColabFold install in **`slurm/run_colabfold_ludmila_rbp.sh`** (and align **`openmpi`** / **`cudasdk`** stacks if RCIC docs require it).

- **Go to project root** (required for **`sbatch`**: Slurm runs the script from a spool copy, so **`SLURM_SUBMIT_DIR`** must be this directory):

  ```bash
  cd /pub/ritwik1/manuscripts/OCPT-Ludmila_PhageTherapy/WTA-RBP_Analysis
  ```
  (Adjust if your clone path differs; always **`cd`** here before **`sbatch slurm/...`**.)

- **Sync only this folder from your laptop** (not the whole `LudmilaCharacterization` tree):

  ```bash
  rsync -avz --progress ./LudmilaCharacterization/analyses/Phage_WTA/ \
    ritwik1@hpc3.rcic.uci.edu:/pub/ritwik1/manuscripts/OCPT-Ludmila_PhageTherapy/WTA-RBP_Analysis/
  ```

  Run from the parent of `LudmilaCharacterization/` (adjust paths if your layout differs).

- Run **`bash bin/sanity_check_paths.sh`** once to create required **`logs/`**, **`data/*`**, **`analysis/`**, and other empty directories used by Slurm and Python scripts.

- **Before first real run:** (1) Confirm **`conda activate phage_wta_env`** works after your **`module load`** / **`conda.sh`** edits in **`config/env_setup.sh`**. (2) Replace **MTESTSEQ** in **`config/wta_reference_genes.fasta`**. (3) Fill in **CUDA/cudnn** **`module load`** lines in **`slurm/run_colabfold_ludmila_rbp.sh`**. (4) Confirm **`colabfold_batch`**, **`prokka`**, and **HMMER** tools are on **PATH** in that env. (5) Test with **`export DRY_RUN=1`** and a **single-array** job; inspect **`logs/`**, then **`unset DRY_RUN`**.

- **Preflight (local):** **`python3 bin/preflight_check.py`** — strict by default (fails on **MTESTSEQ**); use **`--allow-wta-placeholder`** only for smoke tests. **Smoke (no Slurm):** **`bash bin/run_smoke_test.sh`**.
- **Preflight (HPC3 login node):** default **`python3`** may be an old user install (e.g. 3.6). Use **`bash bin/hpc3_login_preflight.sh --allow-wta-placeholder`** (loads an **`anaconda/...`** module, then runs preflight with **`PHAGE_WTA_SKIP_CONDA=1`**), or run **`module load anaconda/2024.06`** yourself before **`python3`**.
- **Local troubleshoot gate (recommended):** **`python3 bin/validate_troubleshoot_complete.py`** — **`bash -n`** on all shell scripts, **`ast.parse`** on Python, verifies Slurm-invoked **`bin/*.py`** exist, runs **`sanity_check_paths`**, **`prepare_genome_lists`**, lenient preflight, and **`build_wta_feature_matrix --allow-empty-raw`**. For an HPC3 command block: **`python3 bin/validate_troubleshoot_complete.py --print-hpc3-commands`**. Customize the printed **`ssh`** line with **`PHAGE_WTA_HPC3_USER`**, **`PHAGE_WTA_HPC3_HOST`**, **`PHAGE_WTA_HPC3_PROJECT_DIR`** (defaults: **`ritwik1`**, **`hpc3.rcic.uci.edu`**, **`/pub/ritwik1/manuscripts/OCPT-Ludmila_PhageTherapy/WTA-RBP_Analysis`**).

- **Preflight (Slurm):** **`sbatch slurm/preflight_check.sh`** — use **`sbatch --export=ALL,PREFLIGHT_ALLOW_WTA_PLACEHOLDER=1`** if references are still placeholders. With **`export DRY_RUN=1`**, the job echoes the **`preflight_check.py`** command and exits **0** without running checks.

---

## 1. Prepare genome lists

- Copy or sync assemblies:
  - Phage: `data/phage/*.fa` or `*.fasta`
  - Hosts: `data/hosts/*.fa` or `*.fasta`
- Generate TSV lists for Slurm:

  ```bash
  python bin/prepare_genome_lists.py
  ```

  Outputs: `config/phage_genomes.tsv`, `config/host_genomes.tsv`.

- If you already populated **`data/host_prokka/<id>/`** (e.g. copied Prokka outputs from **`OCPT-Ludmila_Characterization/outputs/annotation`**) but **`config/host_genomes.tsv`** still has a toy row, rebuild the TSV from directory names and **`../../outputs/assemblies_cleaned/<id>/contigs.fasta`** (manuscript layout on HPC3):

  ```bash
  python bin/sync_host_genomes_from_host_prokka.py
  ```

  Use **`--dry-run`** to print the mapping without writing. Then run **`python3 bin/preflight_check.py --allow-wta-placeholder`** and confirm there are no **host_prokka** vs TSV mismatch warnings.

  If you **`rsync`** this folder from a laptop that still has a toy **`config/host_genomes.tsv`**, re-run **`sync_host_genomes_from_host_prokka.py` on HPC3** after sync so the server keeps all **83** host rows.

---

## 2. Annotate genomes and detect RBPs

- Set **`N`** = number of **data** lines in `config/phage_genomes.tsv` (total lines minus header). Set **`M`** the same way for `config/host_genomes.tsv` (e.g. `N=$(($(wc -l < config/phage_genomes.tsv) - 1))`).

- Submit Prokka:

  ```bash
  sbatch --array=1-N slurm/annotate_phage_prokka.sh
  sbatch --array=1-M slurm/annotate_hosts_prokka.sh
  ```

- After phage Prokka finishes, submit RBP detection (same **`N`** as phage genomes):

  ```bash
  sbatch --array=1-N slurm/detect_phage_rbp.sh
  ```

- Monitor jobs:

  ```bash
  squeue -u "$USER"
  ```

  Wait for relevant jobs to leave the queue before steps that depend on their outputs.

---

## 3. Build WTA HMMs and scan hosts

- Replace placeholder sequences in **`config/wta_reference_genes.fasta`** with curated WTA reference proteins.

- Build the HMM library (once):

  ```bash
  sbatch slurm/build_wta_hmmdb.sh
  ```

- After **host** Prokka is complete, scan hosts for WTA genes (**replace `M`** as in step 1):

  ```bash
  sbatch --array=1-M slurm/scan_hosts_for_wta_genes.sh
  ```

- Merge per-host hit tables, then build the feature matrix:

  ```bash
  python bin/scan_hosts_for_wta_genes.py --collect-raw
  python bin/build_wta_feature_matrix.py
  ```

  Output: **`data/wta_features/wta_feature_matrix.tsv`**.

  Merged **`data/wta_features/raw_wta_hits.tsv`** uses columns **`host_id`**, **`gene`**, **`best_evalue`**, **`bitscore`**, **`hit_length`**. If you have old per-host **`by_host/*_wta_hits.tsv`** from a previous layout, delete them and re-run the scan array (or **`--collect-raw`** only after regenerating **`by_host`** files).

---

## 4. Ludmila RBP extraction and AlphaFold / ColabFold

- Extract candidate RBP/tail proteins (adjust `--prokka-dir` if your phage folder name differs):

  ```bash
  python bin/extract_ludmila_rbp_candidates.py \
    --prokka-dir data/phage_prokka/Ludmila \
    --out data/phage_rbp/Ludmila_RBP_candidates.fasta
  ```

- Submit GPU structure prediction (edit CUDA modules and `colabfold_batch` flags in the script first). The repo uses **`#SBATCH --partition=free-gpu`** so you do **not** need a dedicated GPU Slurm account; jobs may be **preempted**. When your lab has a **`*_gpu`** account, switch the script to **`partition=gpu`** and that account.

  ```bash
  sbatch slurm/run_colabfold_ludmila_rbp.sh
  ```

- ColabFold outputs are written under **`data/alphafold/Ludmila_RBPs/`** (PDBs and run metadata, depending on your ColabFold version and options).

---

## 5. WTA–infectivity association

- Edit **`data/phenotypes/ludmila_infectivity.tsv`** so `host_id` matches **`wta_feature_matrix.tsv`** and `infectivity` reflects your real data (see comment lines in that file).

- Run the association job (or the same Python command interactively):

  ```bash
  sbatch slurm/associate_wta_with_infectivity.sh
  ```

- Regression coefficients and ranks are written to **`analysis/wta_infectivity_logreg_coeffs.tsv`**.

---

## 6. Inspect results

- **`data/wta_features/wta_feature_matrix.tsv`** — per-host WTA gene / backbone / glycan flags.
- **`data/phenotypes/ludmila_infectivity.tsv`** — phenotypes used for the merge.
- **`analysis/wta_infectivity_logreg_coeffs.tsv`** — logistic regression feature coefficients.
- **`data/alphafold/Ludmila_RBPs/`** — predicted structures.

Use **Jupyter** on HPC3, **download** these TSVs and PDBs to your laptop, or open a **local notebook** that reads copies of the tables and structure files.

---

## Dependency reminder

- **Stages 1–2:** Prokka; PhageRBPdetection (if used).
- **Stage 3:** HMMER; Python (scan/matrix scripts).
- **Stage 4:** ColabFold; **`free-gpu`** in repo scripts (or **`gpu`** + GPU lab account); CUDA modules if needed.
- **Stage 5:** pandas, scikit-learn.

If a step fails, check the matching **`logs/`** Slurm output and confirm the conda env was active via `env_setup.sh`.

---

## Dry-run tests

Minimal wiring checks (dummy FASTAs only; no real genomes):

1. **`bash bin/sanity_check_paths.sh`** — ensures expected directories exist.

2. **`python bin/prepare_genome_lists.py`** — rescans **`data/phage/`** and **`data/hosts/`** and refreshes **`config/phage_genomes.tsv`** / **`config/host_genomes.tsv`**. Paths are stored **relative to the project root** when files live under that tree (Slurm jobs resolve them with **`PHAGE_WTA_PROJECT_ROOT`**). The repo ships example rows pointing at **`Ludmila_test.fa`** and **`HostA_test.fa`** for smoke tests.

3. **Prokka Slurm dry-run** — **`DRY_RUN` must be exported** so Slurm forwards it to the job (prefix-only assignment is not always propagated). With **`DRY_RUN=1`**, drivers set **`PHAGE_WTA_SKIP_CONDA=1`** so **`conda activate phage_wta_env`** is not required for echo-only wiring tests (create the env before real runs).

   ```bash
   export DRY_RUN=1
   sbatch --array=1-1 slurm/annotate_phage_prokka.sh
   ```

   Optionally the same for hosts:

   ```bash
   export DRY_RUN=1
   sbatch --array=1-1 slurm/annotate_hosts_prokka.sh
   ```

4. Find the new log under **`logs/`** (e.g. `phage_annotate_*_1.out` for array task 1). Confirm the printed **`prokka ...`** line shows the correct **`--outdir`**, **`--prefix`**, and input FASTA path, and that **`Expected output proteins:`** points to **`.../data/phage_prokka/<ID>/<ID>.faa`** (Prokka’s default when **`--prefix`** matches the genome **`ID`**).

5. **Other Slurm drivers with `DRY_RUN=1`** (same **`export DRY_RUN=1`** before **`sbatch`**): **`slurm/preflight_check.sh`**, **`slurm/detect_phage_rbp.sh`**, **`slurm/build_wta_hmmdb.sh`**, **`slurm/scan_hosts_for_wta_genes.sh`** (array), **`slurm/run_colabfold_ludmila_rbp.sh`**, **`slurm/associate_wta_with_infectivity.sh`**. Each echoes the command or key paths and exits **0** without invoking heavy tools. **`slurm/template_cpu_job.sh`** and **`slurm/template_gpu_job.sh`** echo example commands under **`DRY_RUN=1`** as well.

6. Unset dry-run before real runs: **`unset DRY_RUN`** (or log out).

7. **Local smoke (optional):** **`bash bin/run_smoke_test.sh`** — refreshes genome lists, runs lenient preflight, builds an empty WTA matrix from **`host_genomes.tsv`**; does **not** submit Slurm jobs.

---

## Slurm → helper scripts (sanity)

All **CPU** jobs below use **`#SBATCH --account=katrine_lab`** and **`--partition=standard`**. **`run_colabfold_ludmila_rbp.sh`** uses **`free-gpu`** (see stage 4).

| Slurm script | Invokes |
|--------------|---------|
| `annotate_phage_prokka.sh` | `prokka` (not Python) |
| `annotate_hosts_prokka.sh` | `prokka` |
| `detect_phage_rbp.sh` | `PhageRBPdetection` or `phageRBPdetection.py` |
| `build_wta_hmmdb.sh` | `bin/build_wta_hmmdb.py` |
| `scan_hosts_for_wta_genes.sh` | `bin/scan_hosts_for_wta_genes.py --host-id …` |
| `run_colabfold_ludmila_rbp.sh` | `colabfold_batch` |
| `associate_wta_with_infectivity.sh` | `bin/associate_wta_with_infectivity.py` |
| `preflight_check.sh` | `bin/preflight_check.py` |

Templates **`template_cpu_job.sh`** / **`template_gpu_job.sh`** are placeholders only. **Prokka** with **`--outdir DIR`** and **`--prefix PREFIX`** writes **`PREFIX.faa`**, **`PREFIX.tsv`**, **`PREFIX.gff`**, etc.; pipeline paths **`data/*_prokka/$ID/$ID.faa`** match that default.
