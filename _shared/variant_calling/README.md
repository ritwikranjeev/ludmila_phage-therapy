# Variant calling

The two workflows that produced the mutation calls drawn in Figure 5 and
Supplementary Figure 5. They live here rather than under either figure because
both figures read the same calls, via `_shared/strain_maps/data/mutations/`.

As with `pipeline/`, these are deposited exactly as they ran on the cluster and
are a provenance record rather than a locally runnable workflow.

## `in_vitro_mutants/` — Figure 5a, 5c and Supplementary Figure 5

Illumina sequencing of the phage-resistant mutants and their ancestors, called
against the ancestral assembly with breseq. A Snakemake workflow
(`workflow/Snakefile`) with per-rule conda environments in `envs/` — fastqc,
fastp, bbmap, assembly, quast, checkm, prokka and breseq — and the Slurm entry
points in `slurm_steps/`, numbered in run order. `run_breseq_S16.sh` and
`slurm_steps/04_rescue_s16_highmem.sh` are a re-run of one sample that needed
more memory than the array default.

## `clinical_ont/` — Figure 5b

Oxford Nanopore sequencing of the pre- and post-treatment clinical isolate pair.
A separate Snakemake workflow: QC, Flye assembly, medaka polishing, Clair3
variant calling, then the filtering and table-building steps in `scripts/`.
`scripts/filter_highconf_variants.py` and `scripts/make_mutation_table_manuscript.py`
are the steps that produced the reported variants;
`reports/summary.md` and `tables/mutation_table_for_manuscript.md` are their
output.

## What to check when reading Figure 5

The breseq reference for each comparison determines what the mutations are
called *relative to*. For the in vitro mutants it is the matched ancestor; for
the clinical pair it is the pre-treatment assembly. Both are set in the
respective `config/config.yaml`.
