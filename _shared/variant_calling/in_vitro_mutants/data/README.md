# Input data for Figure 5a mutation pipeline

**Short-read only.** Place your paired-end short reads in **raw_fastqs/** and point **metadata/samples.tsv** at them. Paths in the sample sheet are relative to the pipeline root (the folder that contains `config/`, `workflow/`, `metadata/`).

## Directory layout

```
data/
├── README.md           (this file)
└── raw_fastqs/         Paired-end short reads (R1, R2) — put your inputs here
```

(If you later add long reads for hybrid assembly, you can use a `data/long_reads/` folder and the `long_reads` column in the sample sheet; for short-read-only runs, leave that column blank.)

## raw_fastqs/

Put paired-end FASTQ files here. The pipeline expects one pair per sample, e.g.:

- `{sample_id}_R1.fastq.gz` and `{sample_id}_R2.fastq.gz`

Example for the default **metadata/samples.tsv**:

- **Reference (host):**  
  `USA300.Anc_R1.fastq.gz`, `USA300.Anc_R2.fastq.gz`
- **Mutants:**  
  `S1.USA300.1.6_R1.fastq.gz`, `S1.USA300.1.6_R2.fastq.gz`,  
  `S2.USA300.2.6_R1.fastq.gz`, `S2.USA300.2.6_R2.fastq.gz`, etc.

In **metadata/samples.tsv**, set:

- `read1` = `data/raw_fastqs/{sample_id}_R1.fastq.gz`
- `read2` = `data/raw_fastqs/{sample_id}_R2.fastq.gz`

(Replace `{sample_id}` with the value in the `sample_id` column for that row.)

## After adding data

1. Ensure **metadata/samples.tsv** has one row per sample (host references and resistant mutants) with correct `read1`/`read2` (and optionally `long_reads`) paths under `data/`.
2. Run the pipeline from the pipeline root; validation will check that all listed files exist.
