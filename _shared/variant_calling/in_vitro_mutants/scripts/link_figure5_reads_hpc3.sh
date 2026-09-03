#!/usr/bin/env bash
# Populate Figure5 sequencing/raw with symlinks so metadata/samples.tsv paths resolve.
# On HPC3, ancestors live under $PROJECT/data/raw/<strain>/ and mutants under
# ludmila-characterization/3-res-mut/reads/raw/raw_fastqs/.
#
# Usage (login or compute): PROJECT=/pub/.../OCPT-Ludmila_Characterization ./scripts/link_figure5_reads_hpc3.sh
# Optional: MUT_READS_DIR=... to override mutant FASTQ directory.

set -euo pipefail
PROJECT="${PROJECT:-/pub/ritwik1/manuscripts/OCPT-Ludmila_Characterization}"
DEST="$PROJECT/data/Figure5_Phage_resistant_Mutant_Analysis/sequencing/raw"
ANC="$PROJECT/data/raw"
MUT="${MUT_READS_DIR:-/pub/ritwik1/manuscripts/ludmila-characterization/3-res-mut/reads/raw/raw_fastqs}"

mkdir -p "$DEST"

link_pair() {
  local src="$1" dest="$2"
  if [[ ! -e "$src" ]]; then
    echo "ERROR: missing source: $src" >&2
    return 1
  fi
  ln -sf "$src" "$dest"
}

link_anc() {
  local subdir="$1" stem="$2"
  link_pair "$ANC/$subdir/${stem}_R1.fastq.gz" "$DEST/${stem}_R1.fastq.gz"
  link_pair "$ANC/$subdir/${stem}_R2.fastq.gz" "$DEST/${stem}_R2.fastq.gz"
}

echo "Ancestors from $ANC -> $DEST"
link_anc USA300 S6.USA300.Anc
link_anc OCPT02SA01 S12.OCPT02SA01.Anc
link_anc OCPT04SA01 S18.OCPT04SA01.Anc
link_anc SA412 S24.SA412.Anc
link_anc PedCF S30.PedCF.Anc

echo "Mutants from $MUT -> $DEST"
mutants=(
  S1.USA300.1.6 S2.USA300.2.6 S3.USA300.3.6 S4.USA300.4.6 S5.USA300.5.6
  S7.OCPT02SA01.1.6 S8.OCPT02SA01.2.6 S9.OCPT02SA01.3.6 S10.OCPT02SA01.4.6 S11.OCPT02SA01.5.6
  S13.OCPT04SA01.1.6 S14.OCPT04SA01.2.6 S15.OCPT04SA01.3.6 S16.OCPT04SA01.4.6 S17.OCPT04SA01.5.6
  S19.SA412.1.6 S20.SA412.2.6 S21.SA412.3.6 S22.SA412.4.6 S23.SA412.5.6
  S25.PedCF.1.6 S26.PedCF.2.6 S27.PedCF.3.6 S28.PedCF.4.6 S29.PedCF.5.6
)
for m in "${mutants[@]}"; do
  link_pair "$MUT/${m}_R1.fastq.gz" "$DEST/${m}_R1.fastq.gz"
  link_pair "$MUT/${m}_R2.fastq.gz" "$DEST/${m}_R2.fastq.gz"
done

echo "OK: $(ls -1 "$DEST" | wc -l) filenames under $DEST"
