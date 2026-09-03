# Figure 2

Panels 2c, 2d and 2e are produced by code in this directory. Panels 2a and 2b
are not, for the reason given below.

## Panels 2a and 2b

Panel 2a is a photograph of plaque morphology and panel 2b is a transmission
electron micrograph of the OCPT-Ludmila virion. Both are primary images
recorded at the bench and were prepared for publication by cropping and
scale-bar annotation in image-editing software. Neither panel is produced by
code, and neither involves computation or statistical analysis, so no script or
input file is deposited for them. The virion dimensions reported in the text
were measured from the micrograph in panel 2b.

## Panel 2c — genome map

The genome and its annotation were produced on a compute cluster; the scripts
are in `scripts/2c/genome_pipeline/`, numbered in the order they were run.

| Step | Script | Tool |
|------|--------|------|
| 01 | `01_clean_reads.sh` | fastp 0.23.2, adapter detection, quality trim (window 4, mean Q20), minimum length 60, base correction |
| 02 | `02_assembly.sh` | SPAdes 3.15.4 `--metaviral`, 16 threads, 120 GB |
| 03 | `03_checkv.sh` | CheckV, database `checkv-db-v1.5`, completeness and contamination |
| 04 | `04_pharokka_annotation.sh` | Pharokka, database `pharokka_v1.8.0_databases`, produces `pharokka.gbk` / `pharokka.gff` |
| 05-06 | `05_phold_analysis.sh`, `06_phold_compare.sh` | Phold, structure-aware functional re-annotation of the Pharokka calls |
| 07 | `07_phold_plot.sh` | Phold's own plotting routine |

`genome_pipeline/config.sh` holds the paths and database locations those scripts
read. `genome_pipeline/spades_run_record/` is SPAdes' own record of the run —
`run_spades.sh` is the exact reconstructed command line, and the three YAML files
give the input read manifest and the error-corrected dataset, which together pin
the assembly parameters more precisely than the wrapper script alone.

### Which script drew the deposited panel

Two routes to a circular genome map exist, and it matters which is which.

`output/2c/ludmila_phold_circos.png` and `.svg` — the files deposited here and
used in the manuscript — were written by **step 07, Phold's plotting routine**
(the PNG metadata records matplotlib as the writing library). Regenerating them
requires re-running the Phold steps, which need the Phold databases and a GPU or
CPU structure-prediction pass; they are not reproducible from the files in
`data/2c/` alone.

`scripts/2c/circos_genome_R.R` is a **separate, self-contained R implementation**
of the same idea, using circlize and rtracklayer. It reads `Ludmila.fasta`,
`pharokka.gff` and `pharokka.gbk` — all three of which are deposited in
`data/2c/` — and writes `ludmila_circos.pdf`. It therefore runs from this
repository with no external database, but it produces its own rendering rather
than reproducing the deposited Phold figure.

```
Rscript scripts/2c/circos_genome_R.R <project_root> <output_dir>
```

Requires `install.packages("circlize")` and
`BiocManager::install("rtracklayer")`.

## Panels 2d and 2e

Produced by `scripts/2d/osgc.R` and `scripts/2e/moi_curves.R` from the inputs in
`data/2d/` and `data/2e/`. Run with `./run.sh`. `scripts/2e/verify_against_original.R`
checks that the rewritten MOI script reproduces the original analysis it replaced;
it compares 576 values and currently reports a maximum difference of zero.
