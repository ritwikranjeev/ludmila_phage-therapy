# Upstream pipeline

This directory is the cluster-side workflow that produced the sequence data the
figures are built from: read QC, assembly, annotation, pangenome construction,
defence-system and prophage calling, phenotype preparation, and the pyseer
association runs. `scripts/NN_*.sh` and `slurm/NN_*.slurm` are numbered in the
order they were run, and `config/hpc_env.sh` names the module and conda
environment loaded at each step.

**These scripts are a provenance record, not a runnable local workflow.** They
are deposited exactly as they were executed, with the absolute paths, Slurm
directives, module loads and database locations of the system they ran on. They
will not execute anywhere else without a Slurm scheduler, the same module tree
and local copies of the CheckV, Pharokka, DefenseFinder and Panaroo databases.
They were left unedited deliberately: rewriting the paths would make them
neither runnable nor an accurate record of what was done.

What they are useful for is checking exactly which tool, version and parameters
produced each input file. Recorded versions:

| Step | Tool |
|------|------|
| 02 | fastqc 0.11.9, fastp 0.23.2 |
| 03 | SPAdes 3.15.4, Unicycler 0.5.0, samtools 1.15.1 |
| 04 | QUAST 5.2.0 / 5.3.0, CheckM |
| 04c | minimap2 2.28, bwa 0.7.17, gcc 11.2.0 |
| 05 | Prokka, `--genus Staphylococcus --species aureus` |
| 06 | Panaroo, `--clean-mode strict -t 16 -a core` |
| 07 | IQ-TREE |
| 08-10 | DefenseFinder |
| 11-12 | PhiSpy |
| 14-16, 21 | pyseer, `--pres` with `--distances` (fixed effects, 10 MDS components) |

Conda tool versions were not pinned in `envs/environment.yml`; the environment
names are recorded in `config/hpc_env.sh`. The Panaroo run is Slurm job 50271308
(2026-03-28) and its stdout and stderr are in `logs/`.

## Relationship to the figure directories

Some scripts appear both here and in a figure directory. Where that happens the
figure directory holds the copy relevant to that figure and this directory's
duplicate is excluded from version control:

- `analyses/Figure2_Genomic_Characterization/` -> `fig02/scripts/2c/genome_pipeline/`
- `analyses/Figure5_Phage_resistant_Mutant_Analysis/` -> `_shared/variant_calling/in_vitro_mutants/`
- `analyses/MRSA_ONT_comparative/` -> `_shared/variant_calling/clinical_ont/`

`analyses/Phage_WTA/` is kept: it produced the Prokka annotations of the 83 host
genomes used in Figure 7 and is not duplicated elsewhere.
`analyses/Figure3_Stability_Testing/` and `analyses/Figure4_Physiological_Characterization/`
hold the original R scripts that the versions in `supp02/` and `fig02/`, `fig04/`
replaced; `fig02/scripts/2e/verify_against_original.R` checks its rewrite against
`analyses/Figure4_Physiological_Characterization/4b_Ludmila_MOI/scripts/moi-gc.R`.

Superseded `legacy/` copies are excluded from version control.

## Not included

The sraP-focused association analysis that preceded the Figure 7 reanalysis is
deliberately not deposited. It reported a gene-content association that the
analysis in `fig07/` shows to be an artefact of Panaroo cluster splitting and
clonal structure, and depositing it would document a conclusion the manuscript
does not make.
