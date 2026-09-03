# Stk2-linked Panaroo DIAMOND on HPC3

## Prerequisites on the server

1. Project tree at e.g. `/pub/ritwik1/manuscripts/OCPT-Ludmila_Characterization/` (rsync from laptop).
2. `outputs/pangenome/gene_data.csv` **or** `combined_protein_CDS.fasta` from Panaroo.
3. A DIAMOND protein database, e.g. `uniref50.dmnd` (build with `diamond makedb`).

## Submit

```bash
cd /pub/ritwik1/manuscripts/OCPT-Ludmila_Characterization
mkdir -p logs
export DIAMOND_DB=/path/on/hpc/uniref50.dmnd
# optional: export GENE_DATA=/custom/path/gene_data.csv
sbatch slurm/stk2_linked_panaroo_diamond.slurm
```

Edit `#SBATCH --partition=` if `free` is wrong for your account.

## Resume (parse only)

If DIAMOND already wrote `results/reports/stk2_linked_diamond_raw.tsv`:

```bash
python3 scripts/annotate_stk2_linked_panaroo.py \
  --raw-hits-tsv results/reports/stk2_linked_diamond_raw.tsv \
  --diamond-db /dummy  # not used when --raw-hits-tsv set; still need gene_data for sequences
```

(Use a real `--gene-data-csv` path; `--diamond-db` is ignored when raw hits are passed—fix script if needed.)

Actually current code: if raw_hits_tsv is set, diamond_db is not checked. Good.

## Logs

- `logs/stk2_linked_diamond_<jobid>.out`
- `logs/stk2_linked_diamond_<jobid>.err`
