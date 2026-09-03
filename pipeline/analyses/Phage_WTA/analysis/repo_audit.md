# Phage_WTA repository audit

Generated during scaffold hardening for HPC3 / Slurm (`katrine_lab`). Paths are relative to `analyses/Phage_WTA/`.

## HPC3 troubleshooting session (Apr 2026)

- **Rsync destination** for this project: **`/pub/ritwik1/manuscripts/OCPT-Ludmila_PhageTherapy/WTA-RBP_Analysis/`** (see **`.cursor/rules/hpc3-rsync.mdc`** and **`validate_troubleshoot_complete.py`** defaults).
- **Login-node `python3`:** was **3.6** from a user **`mamba`** install on **`PATH`** before RCIC **`anaconda`** — **`preflight_check.py`** could not parse. **Mitigation:** **`bash bin/hpc3_login_preflight.sh`** (loads **`anaconda/2024.06`** or similar, then runs preflight with **`PHAGE_WTA_SKIP_CONDA=1`**), or **`module load anaconda/2024.06`** before **`python3`**. **`preflight_check.py`** also errors with a clear message if the interpreter is older than **3.9**.
- **`conda env list`** on the cluster did not show **`phage_wta_env`** — create it before Slurm jobs that source **`env_setup.sh`** without **`SKIP_CONDA`**.
- **Slurm spool / `$0`:** first **`sbatch`** attempt failed because **`source "$(dirname "$0")/../config/env_setup.sh"`** resolved to the spool directory, not the repo. **Fix:** all **`slurm/*.sh`** drivers now **`cd "${SLURM_SUBMIT_DIR}"`** (fallback: **`BASH_SOURCE`**) before **`source "${PWD}/config/env_setup.sh"`**, and array jobs source **`_array_tsv_common.sh`** via **`${PHAGE_WTA_PROJECT_ROOT}/slurm/...`**.
- **`DRY_RUN=1` + conda:** Slurm jobs failed with **`EnvironmentNameNotFound: phage_wta_env`** until **`DRY_RUN=1`** triggered **`PHAGE_WTA_SKIP_CONDA=1`** in the drivers so echo-only jobs do not require the env. **Create `phage_wta_env` before any real (non-DRY_RUN) Slurm run.**
- **Preflight (after module load):** HMMER/BLAST on **`PATH`** from user **`mamba`**; **Prokka** and **`colabfold_batch`** still missing until **`phage_wta_env`** (or equivalent) is installed.

## RCIC module inventory (user `module avail`, 2026)

- **Anaconda:** `anaconda/2020.07` … `anaconda/2025.12` under **`LANGUAGES`** — matches **`env_setup.sh`** try-order (override with **`PHAGE_WTA_CONDA_MODULE`**).
- **HMMER / BLAST:** **`BIOTOOLS`** — **`hmmer/3.3`**, **`ncbi-blast/2.10.0`**, **`ncbi-blast/2.13.0`** — optional alternative to installing HMMER/BLAST only inside conda; watch **`PATH`** vs **`phage_wta_env`**.
- **CUDA:** **`TOOLS`** — e.g. **`cuda/11.7.1`**, **`cuda/12.2.0`**, **`cuda/13.0.1`** for ColabFold alignment.
- **Prokka:** not listed as a dedicated module in that snapshot — use **Bioconda** in **`phage_wta_env`**.

## In-place hardening (latest pass)

- **`slurm/preflight_check.sh`**: **`DRY_RUN=1`** echoes the **`preflight_check.py`** invocation and exits **0**.
- **`sanity_check_paths.sh`**: creates **`data/wta_hmms/`** alongside other expected dirs.
- **`raw_wta_hits.tsv` / `by_host/*_wta_hits.tsv`**: columns are **`host_id`**, **`gene`**, **`best_evalue`**, **`bitscore`**, **`hit_length`** (regenerate **`by_host`** if upgrading from the older column order).
- **`config/env_setup.sh`**: clearer stderr hints when **`conda activate phage_wta_env`** fails.
- **`validate_troubleshoot_complete.py`**: **`--print-hpc3-commands`** respects **`PHAGE_WTA_HPC3_USER`**, **`PHAGE_WTA_HPC3_HOST`**, **`PHAGE_WTA_HPC3_PROJECT_DIR`**.

## What looks ready

- **Slurm drivers** (`slurm/*.sh` except templates) all use `#SBATCH --account=katrine_lab`, `set -euo pipefail`, `source config/env_setup.sh`, `cd "${PHAGE_WTA_PROJECT_ROOT}"`, and **`mkdir -p logs`** before writing log files.
- **Genome TSV paths**: `prepare_genome_lists.py` writes paths relative to the project root; **annotate_*** / **detect_phage_rbp** resolve non-absolute paths with **`PHAGE_WTA_PROJECT_ROOT`**.
- **Prokka outputs**: **`--prefix $ID`** + **`--outdir .../$ID`** ⇒ **`$ID.faa`** / **`.tsv`** / **`.gff`** — matches scripts expecting **`data/*_prokka/$ID/$ID.faa`**.
- **Array jobs** share **`slurm/_array_tsv_common.sh`**: validates **`SLURM_ARRAY_TASK_ID`** is in **1..N** for **N** data rows (header excluded).
- **`DRY_RUN=1`**: supported on **preflight_check**, **annotate_phage**, **annotate_hosts**, **detect_phage_rbp**, **build_wta_hmmdb**, **scan_hosts_for_wta_genes**, **run_colabfold_ludmila_rbp**, **associate_wta_with_infectivity**, and on **`template_cpu_job.sh`** / **`template_gpu_job.sh`** (export before **`sbatch`** so the variable reaches the compute node).
- **`config/env_setup.sh`**: sets **`PHAGE_WTA_PROJECT_ROOT`** / **`PHAGE_WTA_DATA_DIR`**; **`PHAGE_WTA_SKIP_CONDA=1`** skips module/conda (used by **`sanity_check_paths.sh`**, **`print_pipeline_status.sh`**, **`run_smoke_test.sh`**).
- **Python layout**: scripts use **`project_root()`** = parent of **`bin/`**; outputs use **`mkdir(parents=True)`** where needed.
- **Preflight**: **`bin/preflight_check.py`** + **`slurm/preflight_check.sh`** check config files, genome TSVs, key dirs, optional tools on **PATH**, and **MTESTSEQ** placeholders (override with **`--allow-wta-placeholder`** or Slurm **`PREFLIGHT_ALLOW_WTA_PLACEHOLDER=1`**).
- **Smoke harness**: **`bin/run_smoke_test.sh`** runs directory setup, genome list refresh, preflight (lenient), and **`build_wta_feature_matrix.py --allow-empty-raw`**.

## What is still stub / TODO

- **PhageRBPdetection** CLI in **`slurm/detect_phage_rbp.sh`** (placeholder two-arg call; tool-specific flags unknown).
- **ColabFold**: **`colabfold_batch`** flags (monomer vs multimer, recycle, etc.) and **GPU/CUDA module** names are comments only.
- **WTA biology**: **`build_wta_feature_matrix.py`** backbone/glycan rules are heuristic placeholders.
- **`extract_ludmila_rbp_candidates.py`**: heuristic RBP picking; Prokka file naming may need adjustment per run.
- **Templates** **`template_cpu_job.sh`** / **`template_gpu_job.sh`**: copy-paste starters only; **`DRY_RUN=1`** echoes example commands, otherwise replace the trailing placeholder echo with real work.

## Likely to fail immediately on HPC3 if not fixed

- **`conda activate phage_wta_env`** when the env does not exist → **`env_setup.sh`** exits **1** (Slurm job fails fast).
- **`module load anaconda/2022.05`** may be wrong or absent → conda not on **PATH** → same failure unless **`conda.sh`** is found under **`$HOME`** or **`PHAGE_WTA_CONDA_MODULE`** is set correctly.
- **`preflight_check.py`** (strict mode) **fails** while **`config/wta_reference_genes.fasta`** still contains **MTESTSEQ** (intentional until you paste real references).
- **RBP / ColabFold / HMMER / Prokka** jobs fail if binaries are not on **PATH** inside the activated env (preflight **warns** but does not fail on missing optional tools except your interpretation of “critical”).
- **`associate_wta_with_infectivity.py`** needs **pandas** + **scikit-learn**; **`ImportError`** exits **1** with a clear message.
- **`detect_phage_rbp.sh`**: will fail if **Prokka** outputs are missing or tool not installed (unless **DRY_RUN=1**).

## Assumptions you must fill manually

- Create **`phage_wta_env`** and install: **Prokka**, **HMMER**, **pandas**, **scikit-learn**, **ColabFold**, **blastp** (if used), **PhageRBPdetection** (if used).
- Set **HPC3 module** names for **anaconda** and **CUDA/cudnn** in **`env_setup.sh`** / **`run_colabfold_ludmila_rbp.sh`**.
- Replace **WTA reference** sequences; align **`host_id`** across **`wta_feature_matrix.tsv`** and **`ludmila_infectivity.tsv`**.
- Confirm **PhageRBPdetection** and **colabfold_batch** invocation and flags.

---

## Current status after automated troubleshooting

*(Updated by the last hardening pass in-repo.)*

- **Code health:** All **`slurm/*.sh`** were **`bash -n`** syntax-checked; all **`bin/*.py`** were checked with **`python3 -m py_compile`** in the maintenance pass (re-run after edits if you change shell/Python).
- **Local (no conda, no Prokka/HMMER):** **`bash bin/run_smoke_test.sh`** uses **`PHAGE_WTA_SKIP_CONDA=1`** and should complete: **`sanity_check_paths`**, **`prepare_genome_lists`**, **`preflight_check.py --allow-wta-placeholder`**, **`build_wta_feature_matrix.py --allow-empty-raw`** (requires **`python3`** + stdlib only for those steps except **numpy** is *not* needed for preflight/prepare/sanity; **build_wta_feature_matrix** uses only stdlib).
- **Strict preflight:** **`python3 bin/preflight_check.py`** is expected to **FAIL** until **MTESTSEQ** is removed from **`config/wta_reference_genes.fasta`** (use **`--allow-wta-placeholder`** only for smoke/tests).
- **HPC3 dry-run:** With a working **`phage_wta_env`**, **`export DRY_RUN=1`** and **`sbatch --array=1-1 slurm/annotate_phage_prokka.sh`** (and host analogue, plus other Slurm drivers listed in **PIPELINE_STEPS**) echo resolved commands/paths and exit **0** without heavy tools. **`slurm/preflight_check.sh`** with **`DRY_RUN=1`** echoes the **`preflight_check.py`** line only.
- **HPC3 preflight Slurm job:** **`sbatch slurm/preflight_check.sh`** runs strict preflight unless you pass **`PREFLIGHT_ALLOW_WTA_PLACEHOLDER=1`** via **`sbatch --export=ALL,PREFLIGHT_ALLOW_WTA_PLACEHOLDER=1`**, or unless you submitted with **`DRY_RUN=1`** (echo-only).

### Troubleshooting continuation (static pass)

- **Slurm → `bin/`:** Every Python helper invoked from **`slurm/*.sh`** exists (**`preflight_check`**, **`build_wta_hmmdb`**, **`scan_hosts_for_wta_genes`**, **`associate_wta_with_infectivity`**). Prokka naming (**`--prefix $ID`** ⇒ **`$ID.faa`**) still matches **`detect_phage_rbp.sh`** and **`extract_ludmila_rbp_candidates.py`**.
- **`run_smoke_test.sh`:** Removed **`|| true`** after **`preflight_check.py --allow-wta-placeholder`** so empty/invalid genome TSVs and other critical preflight errors fail the smoke test instead of being hidden.
- **Slurm Python launcher:** **`preflight_check.sh`**, **`build_wta_hmmdb.sh`**, **`scan_hosts_for_wta_genes.sh`**, and **`associate_wta_with_infectivity.sh`** now call **`python3`** (not **`python`**) so jobs do not break on nodes where only **`python3`** is on **PATH** after **`conda activate`**.

---

## HPC3 preflight fixes

*(Cursor / local repo edits; run the commands below on HPC3 to confirm.)*

1. **`config/env_setup.sh`**: default Anaconda module is now tried in order **`anaconda/2024.06`**, **`anaconda/2022.05`**, **`anaconda/2021.11`** (per [RCIC user-installed / modules documentation](https://rcic.uci.edu/software/user-installed.html)). **`PHAGE_WTA_CONDA_MODULE`** still overrides. If none load, a **WARNING** is printed and the script falls back to **`$HOME` miniforge/miniconda** `conda.sh` if present. **`PHAGE_WTA_SKIP_CONDA=1`** unchanged.
2. **`data/phenotypes/ludmila_infectivity.tsv`**: added **`HostA_test`** row so merges work when the feature matrix is built only from **`config/host_genomes.tsv`** (toy path). Replace with real strains.
3. **`bin/associate_wta_with_infectivity.py`**: minimum merged hosts is configurable via **`--min-hosts`** or **`PHAGE_WTA_ASSOC_MIN_HOSTS`** (default **5**). Use **`2`** for tiny toy association tests.
4. **`slurm/run_colabfold_ludmila_rbp.sh`**: comments added for RCIC **GPU** partitions (**`gpu`**, **`free-gpu`**), **`--gres`** forms, and the fact that **GPU billing accounts may differ** from CPU **`katrine_lab`**—user must confirm with **`sacctmgr`** / PI if jobs are rejected.

**On HPC3 (project root), run:**

```bash
bash bin/sanity_check_paths.sh
python3 bin/preflight_check.py --allow-wta-placeholder
```

If **`conda activate`** fails, create **`phage_wta_env`** on an interactive node (RCIC recommends **`srun`** into **`free`** for conda installs). If **`module load`** fails, set **`export PHAGE_WTA_CONDA_MODULE=...`** to a version from **`module spider anaconda`**.

---

## HPC3 dry-run fixes

**Not executed from Cursor:** this environment cannot **`sbatch`** to HPC3 or read cluster **`logs/`**. After you submit the dry-run batch from **PIPELINE_STEPS** / the user checklist, inspect each **`logs/*.out`** and verify echoed **`prokka`**, **`python .../scan_hosts...`**, **`colabfold_batch`**, etc.

**Repo-side expectation:** relative FASTA paths in **`config/*_genomes.tsv`** are joined with **`PHAGE_WTA_PROJECT_ROOT`** in Slurm; **`DRY_RUN=1`** Prokka blocks echo paths even if FASTA is missing (**NOTE** to stderr). If you see a wrong prefix in logs, check **`cd "${PHAGE_WTA_PROJECT_ROOT}"`** and **`sbatch`** submission **cwd**.

*(Append concrete log findings here after you run dry-jobs on HPC3.)*

---

## Prokka on HPC3

**Pending:** run real **`sbatch --array=1-1`** (phage) and **`1–M`** (hosts) after replacing dummy FASTAs and **`python3 bin/prepare_genome_lists.py`**.

**Current committed TSVs:** **`config/phage_genomes.tsv`** → **1** data row (**`Ludmila_test`**); **`config/host_genomes.tsv`** → **1** data row (**`HostA_test`**). So **`N=1`**, **`M=1`** until you add genomes and regenerate.

*(After jobs: note Prokka versions, failures, and any flag fixes, e.g. **`--metagenome`** for fragmented phage.)*

---

## WTA HMM + scan on HPC3

**Pending:** replace **MTESTSEQ** in **`config/wta_reference_genes.fasta`**, then strict **`python3 bin/preflight_check.py`**, **`sbatch slurm/build_wta_hmmdb.sh`**, **`sbatch --array=1-1 slurm/scan_hosts_for_wta_genes.sh`**, **`python3 bin/build_wta_feature_matrix.py`**.

---

## ColabFold on HPC3

**Pending:** **`run_colabfold_ludmila_rbp.sh`** uses **`free-gpu`** (no **`*_gpu`** account). Still need **ColabFold** install + **`colabfold_batch`** on **PATH**, optional **CUDA** modules, and tolerance for **preemption**. Switch to **`gpu`** + GPU lab account when available.

---

## WTA–infectivity association on HPC3

**Pending:** ensure **`wta_feature_matrix.tsv`** and **`ludmila_infectivity.tsv`** share **`host_id`** values; for fewer than five hosts use **`export PHAGE_WTA_ASSOC_MIN_HOSTS=2`** before **`sbatch slurm/associate_wta_with_infectivity.sh`**, or run **`python .../associate_wta_with_infectivity.py --min-hosts 2`**.

---

## HPC3 run status

- **Minimal Prokka:** Not verified on cluster from Cursor — repo paths and **`katrine_lab`** CPU scripts are consistent; verify with **`DRY_RUN=1`** jobs on HPC3 after **`phage_wta_env`** + **Prokka** exist.
- **Minimal WTA HMM + scan:** Blocked on real WTA FASTA + HMMER in env — strict preflight fails until **MTESTSEQ** is removed.
- **Minimal ColabFold:** Blocked mainly on **ColabFold** install / **`colabfold_batch`**; partition is **`free-gpu`** (no GPU lab account). Preemption possible.
- **Minimal WTA–infectivity association:** Runnable in principle with **`--min-hosts 2`**, matching **`host_id`**, and **pandas/sklearn** in the env.

**Ready for small-scale real execution?** Yes for one phage + one host Prokka once **`phage_wta_env`** exists and **`module load anaconda/...`** succeeds on compute nodes — cluster-side confirmation is your **`validate_troubleshoot_complete.py --print-hpc3-commands`** block plus log inspection.

**Before scaling up:** install/configure **Prokka**, **HMMER**, **ColabFold** (and optional **PhageRBPdetection**, **blastp**), replace **WTA** references, confirm **GPU** Slurm account if using **`run_colabfold_ludmila_rbp.sh`**, and expand **`host_genomes.tsv`** / phenotypes with consistent **`host_id`**s.

---

## Troubleshooting completion (what “done” means)

**In-repo / static (Cursor):** Slurm drivers, **`python3`** invocations, Prokka naming, array validation, **`DRY_RUN`** blocks, and **`bin/`** helpers were re-checked; **`run_smoke_test.sh`** no longer masks preflight failures; **`validate_troubleshoot_complete.py`** was added as the single local gate.

**Local on your machine (required for “local complete”):** From **`analyses/Phage_WTA/`**, run:

```bash
python3 bin/validate_troubleshoot_complete.py
```

Exit **0** ⇒ syntax, wiring, **`sanity_check_paths`**, genome list refresh, lenient preflight, and empty-raw feature matrix all succeeded. (The Cursor agent cannot reliably execute shell/Python against your workspace in this environment, so this command must be run locally.)

**HPC3 (required for “cluster dry-lab complete”):** After **`rsync`** of **`analyses/Phage_WTA/`** to **`/pub/ritwik1/manuscripts/OCPT-Ludmila_PhageTherapy/WTA-RBP_Analysis/`** (or your override path), run the printed checklist:

```bash
python3 bin/validate_troubleshoot_complete.py --print-hpc3-commands
```

Then execute that **`ssh ... bash -lc '...'`** block. Success means: **`preflight_check.py --allow-wta-placeholder`** exits **0** on the login node, and both **`DRY_RUN=1`** Prokka array jobs leave **`logs/phage_annotate_*_1.out`** / **`logs/host_annotate_*_1.out`** containing **`[DRY_RUN]`** and the echoed **prokka** command lines. Full scientific tracks (real Prokka, HMM, ColabFold, association) remain **manual** until tools, references, and GPU account are configured — see **HPC3 run status** above.

**`validate_troubleshoot_complete.py --print-hpc3-commands`:** override SSH user, host, and remote project directory with **`PHAGE_WTA_HPC3_USER`**, **`PHAGE_WTA_HPC3_HOST`**, **`PHAGE_WTA_HPC3_PROJECT_DIR`**.

**Evidence in this session:** Your Cursor terminal showed **SSH to HPC3 working** (**`login-i17`**), but **no saved transcript** here of **`sbatch`** dry-runs against the synced tree — mark HPC3 dry-lab **complete** only after you run the printed block and see those log files.
