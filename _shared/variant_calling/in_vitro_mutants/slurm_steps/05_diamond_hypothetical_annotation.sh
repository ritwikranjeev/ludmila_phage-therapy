#!/bin/bash
#SBATCH -A katrine_lab
#SBATCH -p standard
#SBATCH -t 02:00:00
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -J diamond_hypo
#SBATCH -o logs/05_diamond_hypo_%j.out
#SBATCH -e logs/05_diamond_hypo_%j.err
#
# One-time: download UniProt *Staphylococcus aureus* (taxon 1280) proteins, build DIAMOND DB.
# Then: extract hypothetical proteins from each host-reference .faa and run diamond blastp.

set -euo pipefail
ANALYSIS_DIR="/pub/ritwik1/manuscripts/OCPT-Ludmila_Characterization/analyses/Figure5_Phage_resistant_Mutant_Analysis"
cd "$ANALYSIS_DIR"
mkdir -p logs

module load diamond/2.0.15

RES="${ANALYSIS_DIR}/resources/diamond_s_aureus_uniprot"
FASTA="${RES}/s_aureus_taxon1280_uniprot.fasta"
DB="${RES}/s_aureus_taxon1280"
QUERY_ROOT="${ANALYSIS_DIR}/results/annotation_diamond_hypotheticals"
UNIPROT_URL='https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(taxonomy_id:1280)+AND+(fragment:false)'

mkdir -p "$RES" "$QUERY_ROOT"

echo "=== DIAMOND hypothetical annotation $(date) ==="
echo "Node: $(hostname)"

if [[ ! -s "${DB}.dmnd" ]]; then
  echo "Downloading UniProt FASTA (taxon 1280, non-fragment) ..."
  mkdir -p "$RES"
  curl -fsSL "$UNIPROT_URL" -o "$FASTA"
  echo "Building DIAMOND database (this may take several minutes) ..."
  diamond makedb --in "$FASTA" -d "$DB" --threads "${SLURM_CPUS_PER_TASK:-8}"
else
  echo "Using existing DIAMOND DB: ${DB}.dmnd"
fi

shopt -s nullglob
refs=( "${ANALYSIS_DIR}/results/references/"*.Anc/*.faa )
if [[ ${#refs[@]} -eq 0 ]]; then
  echo "ERROR: No reference .faa under results/references/*/*.faa" >&2
  exit 1
fi

for faa in "${refs[@]}"; do
  base=$(basename "$faa" .faa)
  outdir="${QUERY_ROOT}/${base}"
  mkdir -p "$outdir"
  hypo="${outdir}/${base}_hypothetical.faa"
  tsv="${outdir}/${base}_diamond_best.tsv"

  echo "--- $base ---"
  python3 "${ANALYSIS_DIR}/scripts/extract_hypothetical_proteins.py" "$faa" "$hypo"
  if [[ ! -s "$hypo" ]]; then
    echo "No hypothetical proteins for $base; skipping DIAMOND."
    echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" > "$tsv"
    continue
  fi

  diamond blastp \
    -d "${DB}.dmnd" \
    -q "$hypo" \
    -o "$tsv" \
    --threads "${SLURM_CPUS_PER_TASK:-8}" \
    --max-target-seqs 1 \
    --max-hsps 1 \
    --evalue 1e-5 \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

  {
    echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"
    cat "$tsv"
  } > "${tsv}.with_header"
  mv "${tsv}.with_header" "$tsv"

  n=$(wc -l < "$tsv")
  n=$((n - 1))
  echo "Wrote $tsv ($n hits)"
done

echo "=== Done $(date) ==="
