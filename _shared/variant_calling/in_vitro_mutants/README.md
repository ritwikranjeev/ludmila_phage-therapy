# Figure 5 – Phage-resistant Mutant Analysis

## Panel A: Mutation-calling pipeline (Snakemake, HPC3-ready)

Production pipeline for resistant-mutant mutation calling: reference assembly, annotation, read QC, breseq, mutation tables, and Figure 5a mutation-table visualization.

### Syncing to HPC3

From your **local machine**, sync the **whole project** (recommended):

```bash
# Option 1: from inside LudmilaCharacterization, run the project sync script
cd /path/to/LudmilaCharacterization
./sync_to_hpc3.sh

# Option 2: from the folder that contains LudmilaCharacterization (e.g. cursor/)
cd /path/to/cursor
rsync -avz --progress ./LudmilaCharacterization/ ritwik1@hpc3.rcic.uci.edu:/pub/ritwik1/manuscripts/OCPT-Ludmila_Characterization/
```

Use the path that contains the **entire** `LudmilaCharacterization` folder (e.g. `cursor/` or your repo root). Do not run rsync from inside `LudmilaCharacterization` with `./Figure5_...` — that path does not exist there (Figure5 is under `analyses/`).

After rsync, on HPC3 set the project root so paths resolve correctly. Either:

- Run Snakemake from the **pipeline root** (the folder that contains `config/`, `workflow/`, `metadata/`), and leave `project_dir: "."` in `config/config.yaml`, or  
- Override when running:  
  `snakemake -s workflow/Snakefile --config project_dir=$PWD ...`

### HPC3 setup (one-time)

SSH to HPC3 and load modules, then create conda envs for tools not provided as modules:

```bash
# Load HPC3 modules (adjust versions if needed)
module load python/3.10.2 miniconda3/24.9.2 mamba/24.3.0
module load fastp/0.23.2 fastqc/0.11.9
module load SPAdes/3.15.4 unicycler/0.5.0
module load samtools/1.15.1 minimap2/2.28 bowtie2/2.5.1 pigz/2.6
module load R/4.4.2

# Conda envs for breseq, Prokka, QUAST (no HPC3 module). CheckM is optional (off by default).
cd /path/to/Figure5_Phage_resistant_Mutant_Analysis
mamba env create -f envs/breseq.yaml
mamba env create -f envs/quast.yaml
# If Prokka fails (parallel dependency conflict), use flexible channel priority:
mamba env create -f envs/prokka.yaml --channel-priority flexible
# CheckM is optional: config has run_checkm: false by default. If the CheckM env fails to solve (pysam/libdeflate conflicts), leave run_checkm: false and the pipeline skips CheckM. If you want CheckM, try:
#   mamba env create -f envs/checkm.yaml --channel-priority flexible
#   then set run_checkm: true and run with mamba activate checkm (Snakemake does not use conda for this rule).
# BBMap: use module load bbmap if available; otherwise mamba env create -f envs/bbmap.yaml for PhiX removal (optional).
```

Alternatively, let Snakemake create envs at run time: use `--use-conda --conda-prefix envs/conda` when you run Snakemake; it will create envs from the yaml files. If Prokka creation fails, run the `mamba env create -f envs/prokka.yaml --channel-priority flexible` command above first. **CheckM is not created by Snakemake** (rule has no conda); leave `run_checkm: false` unless you install CheckM yourself and activate it before running.

QUAST is provided via the `envs/quast.yaml` conda env (create it with `mamba env create -f envs/quast.yaml`).

**Using HPC3 modules for fastp / fastqc:** If conda env creation for fastp fails (e.g. `LibMambaUnsatisfiableError` for libdeflate), the pipeline can use HPC3 modules instead. In `config/config.yaml`, the `modules:` section sets **fastp**, **fastqc**, and optionally **unicycler**. The **gcc** module (e.g. `gcc/11.2.0`) is loaded with fastp so jobs get a libstdc++ that provides `GLIBCXX_3.4.26`. If SLURM jobs still fail, set **modules.gcc_lib_path** to your GCC 11 lib64 path (e.g. `/opt/apps/gcc/11.2.0/lib64`) so the rule prepends it to `LD_LIBRARY_PATH` before running fastp.

### Config (edit for your allocation)

Edit `config/config.yaml`:

- **project_dir:** Leave as `"."` if you run from the pipeline root; or set to your full HPC3 path (e.g. `/pub/ritwik1/.../Figure5_Phage_resistant_Mutant_Analysis`).
- **cluster.account**, **cluster.partition:** Set to your SLURM account and partition (e.g. `katrine_lab`, `standard`).
- **raw_reads_dir**, **samples_tsv:** Paths relative to `project_dir` (or override via `--config`).
- **exclude_samples:** Optional list of `sample_id` to skip for breseq and downstream (e.g. `["S16.OCPT04SA01.4.6"]` if one sample repeatedly OOMs). Pipeline runs for the remaining mutants.

### Sample sheet

- **metadata/samples.tsv** must have columns: `sample_id`, `sample_type`, `host_id`, `reference_group`, `read1`, `read2` (and optionally `long_reads`, `replicate`, `notes`).
- Set `sample_type` to `host_reference` for reference genomes and `resistant_mutant` for mutants. Each mutant’s `reference_group` must match a `sample_id` of a `host_reference` row.
- Paths in `read1`/`read2` are relative to `project_dir` unless absolute.

### Validate sample sheet and paths

From the **pipeline root** (directory containing `config/`, `workflow/`, `metadata/`):

```bash
python3 scripts/validate_samples.py metadata/samples.tsv
```

(Use `--project-dir /path/to/pipeline/root` if running from elsewhere.)

### Snakemake (if not already available)

If `snakemake` is not found after loading modules, use a dedicated env:

```bash
# From pipeline root, with miniconda/mamba loaded:
mamba env create -f envs/snakemake.yaml
mamba activate snakemake
```

Then run the dry-run and pipeline commands below. (If you installed Snakemake in base, run `conda activate base` so `snakemake` is on PATH.)

### Dry-run and local test

```bash
# Dry-run (see planned jobs)
snakemake -s workflow/Snakefile -n

# Local run (single thread, no SLURM)
snakemake -s workflow/Snakefile --cores 4
```

### SLURM / HPC3 execution (Snakemake 8+)

Snakemake 8+ removed `--cluster` and `--cluster-config`. Use the SLURM executor plugin instead.

**1. Install the SLURM executor plugin** (in your snakemake env):

```bash
mamba activate snakemake
pip install snakemake-executor-plugin-slurm
```

**2. Satisfy Snakemake’s conda version (≥ 24.7.1)**  
Snakemake requires conda ≥ 24.7.1 for `--use-conda`. If your system/module conda is older (e.g. 24.3.0) and you cannot update it, install a **separate** Miniconda in your home and use it only when running Snakemake:

```bash
# One-time: install Miniconda3 in your home (does not touch system conda)
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3_snakemake
rm Miniconda3-latest-Linux-x86_64.sh
```

Then **whenever you run Snakemake**, put this conda first in PATH so Snakemake uses it:

```bash
export PATH=$HOME/miniconda3_snakemake/bin:$PATH
cd /pub/ritwik1/manuscripts/OCPT-Ludmila_Characterization/analyses/Figure5_Phage_resistant_Mutant_Analysis
mamba activate snakemake   # still use your existing snakemake env for the snakemake command
snakemake -s workflow/Snakefile --profile profile/slurm --config project_dir=$PWD
```

**3. Edit `profile/slurm/config.yaml`** and set `slurm_account` and `slurm_partition` to your allocation (same as `config/config.yaml` cluster section).

**4. Run with the profile:**

```bash
mkdir -p logs/slurm
# If using the separate Miniconda above, run this first:
#   export PATH=$HOME/miniconda3_snakemake/bin:$PATH
snakemake -s workflow/Snakefile --profile profile/slurm --use-conda --conda-prefix envs/conda --config project_dir=$PWD
```

**4b. Run via SLURM (survives logout):** Submit the driver as a job so you can log out while the workflow runs:

```bash
cd /pub/ritwik1/manuscripts/OCPT-Ludmila_Characterization/analyses/Figure5_Phage_resistant_Mutant_Analysis
sbatch run_snakemake.sh
```

**Tracking your runs**

- **SLURM job queue** – See what’s running or pending:
  ```bash
  squeue -u $USER
  ```
  Columns: JOBID, PARTITION, NAME, USER, ST, TIME, NODES, NODELIST(REASON). `snakemake_fig5` is the driver; other names are rule names (e.g. breseq, ref_prokka).

- **Driver job output** (when you used `sbatch run_snakemake.sh`):
  ```bash
  cd /pub/ritwik1/manuscripts/OCPT-Ludmila_Characterization/analyses/Figure5_Phage_resistant_Mutant_Analysis
  tail -f logs/snakemake_<JOBID>.out   # live progress (steps X of 58 done, jobs submitted)
  tail logs/snakemake_<JOBID>.err      # errors from the driver
  ```
  Replace `<JOBID>` with the ID from `sbatch` (e.g. `49947125`). Use `ls -t logs/snakemake_*.out` to find the latest.

- **Recent job history** (finished/cancelled jobs, exit codes, runtime):
  ```bash
  sacct -u $USER --starttime today
  ```

- **Per-rule SLURM logs** (when a rule fails, check the log path in the error message):
  ```bash
  ls .snakemake/slurm_logs/                    # by rule (e.g. rule_breseq, rule_ref_prokka)
  cat .snakemake/slurm_logs/rule_breseq/S25.PedCF.1.6/<jobid>.log   # example
  ```

- **Snakemake run log** (one per run, full trace):
  ```bash
  ls -t .snakemake/log/*.snakemake.log | head -1
  ```

The profile sets `executor: slurm`, **max 4 concurrent jobs** (`jobs: 4`), `use-conda: true`, and default SLURM resources. **fastp** and **ref_assembly** get extra time and CPUs via `slurm_extra` (4h/16 CPUs and 8h/16 CPUs). To override job count: `snakemake ... --profile profile/slurm --jobs 10 ...`.

**If a rule fails with `message: None`:** check the SLURM log for that job (path is printed in the error). For example:
`cat .snakemake/slurm_logs/rule_fastp/<sample_id>/<jobid>.log`
to see the real error (e.g. timeout, GLIBCXX, or command failure).

**Cancelling a run**
- In the terminal where Snakemake is running: press **Ctrl+C** once (or twice if needed). Snakemake will shut down and stop submitting new jobs.
- To cancel your SLURM jobs as well: `scancel -u $USER` (all your jobs) or `scancel <jobid>` for a single job. List jobs with `squeue -u $USER`.

**Diagnosing TIMEOUTs**  
If you see `SLURM status is: 'TIMEOUT'` and job resources show `runtime=4, cpus=1`, the executor was requesting only ~4 minutes and 1 CPU. **Time and CPUs are now set in the Snakefile** (e.g. runtime=1440 for breseq). Re-run after pulling the latest workflow.

**Diagnosing FAILED (non-TIMEOUT) jobs**  
If you see `SLURM status is: 'FAILED'`, the job exited with an error (e.g. breseq crash, OOM). Check the rule log path in the message, e.g.  
`cat .snakemake/slurm_logs/rule_breseq/<sample_id>/<jobid>.log`  
Fix the cause then re-run; Snakemake will only re-run the failed job and downstream steps.

**Compute-node I/O error (e.g. one sample always FAILED)**  
If the rule log shows `OSError: [Errno 5] Input/output error: '/data/homezvol1/...'` (Snakemake trying to create its source cache on the compute node’s home), the node’s home filesystem had an I/O error. The pipeline now sets `XDG_CACHE_HOME` to the project’s `.snakemake/xdg_cache` in `run_snakemake.sh` so the cache lives on project storage instead of home. Re-run with that change; if you run Snakemake directly (not via `sbatch run_snakemake.sh`), set `export XDG_CACHE_HOME=$PWD/.snakemake/xdg_cache` before running.

**Partition memory limit (standard)**  
On HPC3, `standard` has `MaxMemPerCPU=6144` MB, so with 16 CPUs the per-job maximum is ~98 GB. Requesting more (e.g. 120 GB) may be capped or fail. For breseq you can set `mem_mb: 96000` in config/resources or the Snakefile to stay under the limit.

**Time limits:** You must request a limit; clusters do not allow unlimited time. The rule-level `slurm_extra` uses 4h for fastp and 24h for assembly. Using a high limit (e.g. 24h) is not risky—jobs still stop when they finish or when the limit is reached.

### Pipeline outputs

- **results/validation/passed.flag** – Sample sheet validation.
- **results/qc/cleaned/** – Cleaned reads per sample (fastp).
- **results/assemblies/{ref_id}.fasta** – Reference assemblies; **references/{ref_id}/** – Prokka annotation (e.g. `.gbk` for breseq).
- **results/assemblies/{ref_id}_qc_compiled.txt** – Compiled QC (QUAST, CheckM if run_checkm, seqkit).
- **results/breseq/{sample_id}/** – Full breseq output per mutant.
- **results/tables/mutations/{sample_id}.tsv** – Per-mutant mutation TSV (gdtools).
- **results/tables/figure5a_plot_table.tsv** – Normalized table for plotting.
- **results/figures/figure5a_mutation_table.png**, **.pdf** – Figure 5a mutation table plot.

---

## Panel C: Kirby Bauer

- **Script:** `scripts/kirby.R`
- **Input:** `data/Figure5_Phage_resistant_Mutant_Analysis/kirby_bauers/5c_KirbyBauer_Input.tsv`
- **Output:** `results/5c_KirbyBauer.png`

Zone diameter (mm) by strain and antibiotic, comparing Ancestor vs Resistant clones. Box-and-whisker plot with overlaid points; Wilcoxon test per Strain × Antibiotic, BH correction; bracket and significance stars for significant differences. Strain order: USA300, SA412, OCPT02SA01, OCPT04SA01, PedCF. Figure dimensions: 22 × 5.5 in (one row of 5 strain panels).

Manuscript figure: run `Rscript manuscript/figures/build_fig5.R` from project root; output `manuscript/figures/Figure5/Fig5.png`.
