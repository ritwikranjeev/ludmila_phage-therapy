#!/usr/bin/env python3
"""
Compute assembly metrics from FASTA (and optionally QUAST report).
Output: TSV with total_length, contig_count, N50, GC_percent, largest_contig, etc.
"""
import argparse
import sys
from pathlib import Path


def parse_fasta(path):
    """Yield (seq_id, seq) for each sequence."""
    name = None
    buf = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if name is not None:
        yield name, "".join(buf)


def compute_metrics(fasta_path, quast_report_path=None):
    lengths = []
    gc_count = 0
    atgc = 0
    for _id, seq in parse_fasta(fasta_path):
            L = len(seq)
            lengths.append(L)
            s = seq.upper()
            gc_count += s.count("G") + s.count("C")
            atgc += s.count("A") + s.count("T") + s.count("G") + s.count("C")
    lengths.sort(reverse=True)
    total = sum(lengths)
    n50 = None
    cum = 0
    for L in lengths:
        cum += L
        if cum >= total / 2:
            n50 = L
            break
    gc_pct = (100.0 * gc_count / atgc) if atgc else 0
    largest = lengths[0] if lengths else 0
    metrics = {
        "total_length": total,
        "contig_count": len(lengths),
        "N50": n50 or 0,
        "GC_percent": round(gc_pct, 2),
        "largest_contig": largest,
    }
    if quast_report_path and Path(quast_report_path).exists():
        with open(quast_report_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 2:
                    key = parts[0].strip()
                    val = parts[1].strip()
                    if key and val and key not in metrics:
                        try:
                            metrics[key] = float(val) if "." in val else int(val)
                        except ValueError:
                            metrics[key] = val
    return metrics


def main():
    p = argparse.ArgumentParser(description="Assembly metrics from FASTA")
    p.add_argument("fasta", help="Assembly FASTA")
    p.add_argument("-o", "--output", required=True, help="Output TSV path")
    p.add_argument("--quast-report", default=None, help="Optional QUAST report.txt")
    args = p.parse_args()

    metrics = compute_metrics(args.fasta, args.quast_report)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        f.write("\t".join(metrics.keys()) + "\n")
        f.write("\t".join(str(v) for v in metrics.values()) + "\n")
    print(f"Wrote {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
