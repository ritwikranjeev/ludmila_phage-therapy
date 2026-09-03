#!/usr/bin/env python3
"""
Per-sample contamination diagnostics:
  1) Map reads to USA300 chromosome (CP000046.1) — flagstat + stats.
  2) Map reads to a 2-sequence panel (Staph chr + E. coli K-12 chr) — count primary hits per contig.
  3) Map reads to assembly contigs — samtools coverage + merge per-contig GC%.

Reads paths come from config/samples.tsv (long or illumina).
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import ssl
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.request import urlopen


def read_url(url: str, timeout: int = 180) -> str:
    try:
        return urlopen(url, timeout=timeout, context=ssl.create_default_context()).read().decode()
    except ssl.SSLError:
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        return urlopen(url, timeout=timeout, context=ctx).read().decode()


def run(cmd: List[str], cwd: Optional[Path] = None) -> None:
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write(r.stdout or "")
        sys.stderr.write(r.stderr or "")
        sys.exit(f"ERROR ({r.returncode}): {' '.join(cmd)}")


def run_out(cmd: List[str], cwd: Optional[Path] = None) -> str:
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write(r.stderr or "")
        sys.exit(f"ERROR ({r.returncode}): {' '.join(cmd)}")
    return r.stdout


def load_sample_row(samples_tsv: Path, sample: str) -> Tuple[str, Optional[Path], Optional[Path], Optional[Path]]:
    with samples_tsv.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("sample") != sample:
                continue
            plat = (row.get("platform") or "").strip()
            r1 = row.get("r1", "")
            r2 = row.get("r2", "")
            lr = row.get("longread", "")
            p1 = Path(r1) if r1 and r1 != "NA" else None
            p2 = Path(r2) if r2 and r2 != "NA" else None
            plr = Path(lr) if lr and lr != "NA" else None
            return plat, p1, p2, plr
    sys.exit(f"ERROR: sample {sample!r} not in {samples_tsv}")


def ensure_fasta(url: str, dest: Path) -> None:
    if dest.is_file() and dest.stat().st_size > 1000:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {dest.name} …", file=sys.stderr)
    dest.write_text(read_url(url, timeout=180))


def _first_fasta_sequence(fa_text: str) -> str:
    """Concatenate sequence lines for the first record in a FASTA string."""
    seq_chunks: List[str] = []
    seen_header = False
    for line in fa_text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if seen_header:
                break
            seen_header = True
            continue
        if seen_header:
            seq_chunks.append(line)
    return "".join(seq_chunks)


def ensure_panel(panel_path: Path) -> None:
    if panel_path.is_file() and panel_path.stat().st_size > 10000:
        return
    panel_path.parent.mkdir(parents=True, exist_ok=True)
    staph_u = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        "db=nucleotide&id=CP000046.1&rettype=fasta&retmode=text"
    )
    ecoli_u = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        "db=nucleotide&id=NC_000913.3&rettype=fasta&retmode=text"
    )
    s = read_url(staph_u, timeout=180)
    e = read_url(ecoli_u, timeout=180)

    def wrap(name: str, seq: str) -> str:
        lines = [f">{name}"]
        for j in range(0, len(seq), 60):
            lines.append(seq[j : j + 60])
        return "\n".join(lines) + "\n"

    a = wrap("Staphylococcus_aureus_CP000046", _first_fasta_sequence(s))
    b = wrap("Escherichia_coli_NC000913", _first_fasta_sequence(e))
    panel_path.write_text(a + b)
    print(f"Wrote panel ({panel_path.stat().st_size} bytes)", file=sys.stderr)


def read_fasta_seqs(path: Path) -> Dict[str, str]:
    out: Dict[str, List[str]] = {}
    cur: Optional[str] = None
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                cur = line[1:].split()[0]
                out[cur] = []
            elif cur:
                out[cur].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


def gc_pct(seq: str) -> float:
    if not seq:
        return 0.0
    s = seq.upper()
    g = s.count("G") + s.count("C")
    return 100.0 * g / len(s)


def map_long_reads_to_bam(ref: Path, reads: Path, bam: Path, threads: int) -> None:
    p1 = subprocess.run(
        [
            "minimap2",
            "-ax",
            "map-ont",
            "-t",
            str(threads),
            str(ref),
            str(reads),
        ],
        capture_output=True,
    )
    if p1.returncode != 0:
        sys.stderr.write(p1.stderr.decode())
        sys.exit("ERROR: minimap2 failed")
    sam = bam.with_suffix(".sam")
    sam.write_bytes(p1.stdout)
    run(
        ["samtools", "sort", "-@", str(threads), "-o", str(bam), str(sam)]
    )
    sam.unlink(missing_ok=True)
    run(["samtools", "index", str(bam)])


def map_pe_to_bam(ref: Path, r1: Path, r2: Path, bam: Path, work: Path, threads: int) -> None:
    ref_copy = work / "ref.fa"
    shutil.copy2(ref, ref_copy)
    run(["bwa", "index", str(ref_copy.name)], cwd=work)
    sam = work / "pe.sam"
    p1 = subprocess.run(
        [
            "bwa",
            "mem",
            "-t",
            str(threads),
            str(ref_copy.name),
            str(r1.resolve()),
            str(r2.resolve()),
        ],
        cwd=work,
        capture_output=True,
    )
    if p1.returncode != 0:
        sys.stderr.write(p1.stderr.decode())
        sys.exit("ERROR: bwa mem failed")
    sam.write_bytes(p1.stdout)
    run(["samtools", "sort", "-@", str(threads), "-o", str(bam), str(sam)])
    sam.unlink(missing_ok=True)
    run(["samtools", "index", str(bam)])
    for suf in (".amb", ".ann", ".bwt", ".pac", ".sa"):
        (work / (ref_copy.name + suf)).unlink(missing_ok=True)


def map_contigs_long(contigs: Path, reads: Path, bam: Path, work: Path, threads: int) -> None:
    ctg = work / "contigs.fa"
    shutil.copy2(contigs, ctg)
    map_long_reads_to_bam(ctg, reads, bam, threads)


def map_contigs_pe(contigs: Path, r1: Path, r2: Path, bam: Path, work: Path, threads: int) -> None:
    ctg = work / "contigs.fa"
    shutil.copy2(contigs, ctg)
    run(["bwa", "index", str(ctg.name)], cwd=work)
    sam = work / "ctg_pe.sam"
    p1 = subprocess.run(
        [
            "bwa",
            "mem",
            "-t",
            str(threads),
            str(ctg.name),
            str(r1.resolve()),
            str(r2.resolve()),
        ],
        cwd=work,
        capture_output=True,
    )
    if p1.returncode != 0:
        sys.stderr.write(p1.stderr.decode())
        sys.exit("ERROR: bwa mem contigs failed")
    sam.write_bytes(p1.stdout)
    run(["samtools", "sort", "-@", str(threads), "-o", str(bam), str(sam)])
    sam.unlink(missing_ok=True)
    run(["samtools", "index", str(bam)])
    for suf in (".amb", ".ann", ".bwt", ".pac", ".sa"):
        (work / ("contigs.fa" + suf)).unlink(missing_ok=True)


def write_flagstat(bam: Path, dest: Path) -> None:
    txt = run_out(["samtools", "flagstat", str(bam)])
    dest.write_text(txt)


def write_stats(bam: Path, dest: Path) -> None:
    txt = run_out(["samtools", "stats", str(bam)])
    dest.write_text(txt)


def panel_hit_counts(bam: Path, dest: Path) -> None:
    # primary alignments only, mapped
    p = subprocess.run(
        [
            "samtools",
            "view",
            "-F",
            "260",
            str(bam),
        ],
        capture_output=True,
        text=True,
    )
    if p.returncode != 0:
        sys.exit("ERROR: samtools view panel")
    counts: Dict[str, int] = {}
    for line in p.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        rname = parts[2]
        if rname == "*":
            continue
        counts[rname] = counts.get(rname, 0) + 1
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["panel_contig", "primary_alignments"])
        for k in sorted(counts.keys()):
            w.writerow([k, counts[k]])


def samtools_coverage(bam: Path, tsv: Path) -> None:
    run(["samtools", "coverage", "-o", str(tsv), str(bam)])


def merge_gc_coverage(
    contigs_fa: Path, cov_tsv: Path, merged: Path
) -> None:
    seqs = read_fasta_seqs(contigs_fa)
    hdr: Optional[List[str]] = None
    cov_rows: List[List[str]] = []
    with cov_tsv.open() as fh:
        for line in fh:
            if line.startswith("#"):
                hdr = line[1:].strip().split("\t")
                continue
            if line.strip():
                cov_rows.append(line.rstrip("\n").split("\t"))
    if not hdr:
        hdr = [
            "rname",
            "startpos",
            "endpos",
            "numreads",
            "covbases",
            "coverage",
            "meandepth",
            "meanbaseq",
            "meanmapq",
        ]
    if not cov_rows:
        merged.write_text("no_coverage_data\n")
        return
    with merged.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(hdr + ["contig_length_fasta", "gc_percent"])
        for parts in cov_rows:
            if len(parts) < 7:
                continue
            cid = parts[0]
            slen = len(seqs.get(cid, ""))
            gcp = gc_pct(seqs.get(cid, ""))
            w.writerow(parts + [str(slen), f"{gcp:.2f}"])


def parse_flagstat_mapped_pct(text: str) -> Tuple[Optional[int], Optional[int], Optional[float]]:
    total = None
    mapped = None
    for line in text.splitlines():
        m = re.match(r"^(\d+) \+ (\d+) in total", line)
        if m:
            total = int(m.group(1)) + int(m.group(2))
        m = re.match(r"^(\d+) \+ (\d+) mapped", line)
        if m:
            mapped = int(m.group(1)) + int(m.group(2))
    pct = 100.0 * mapped / total if total and mapped is not None and total > 0 else None
    return total, mapped, pct


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--samples-tsv", type=Path, required=True)
    ap.add_argument("--assembly", type=Path, required=True)
    ap.add_argument("--usa300-ref", type=Path, required=True)
    ap.add_argument("--panel-ref", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    plat, r1, r2, lr = load_sample_row(args.samples_tsv, args.sample)
    ensure_fasta(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        "db=nucleotide&id=CP000046.1&rettype=fasta&retmode=text",
        args.usa300_ref,
    )
    ensure_panel(args.panel_ref)

    if not args.assembly.is_file():
        sys.exit(f"ERROR: missing assembly {args.assembly}")

    out = args.out_dir / args.sample
    out.mkdir(parents=True, exist_ok=True)
    tmp_root = os.environ.get("TMPDIR") or os.environ.get("SLURM_TMPDIR") or None
    work = Path(tempfile.mkdtemp(prefix=f"diag_{args.sample}_", dir=tmp_root))

    try:
        # --- USA300 reference ---
        usa_bam = work / "to_usa300.bam"
        pan_bam = work / "to_panel.bam"
        ctg_bam = work / "to_contigs.bam"

        if lr and lr.is_file():
            map_long_reads_to_bam(args.usa300_ref, lr, usa_bam, args.threads)
            map_long_reads_to_bam(args.panel_ref, lr, pan_bam, args.threads)
            map_contigs_long(args.assembly, lr, ctg_bam, work, args.threads)
        elif r1 and r2 and r1.is_file() and r2.is_file():
            map_pe_to_bam(args.usa300_ref, r1, r2, usa_bam, work, args.threads)
            map_pe_to_bam(args.panel_ref, r1, r2, pan_bam, work, args.threads)
            map_contigs_pe(args.assembly, r1, r2, ctg_bam, work, args.threads)
        else:
            sys.exit("ERROR: no reads for sample")

        write_flagstat(usa_bam, out / "map_to_USA300.flagstat.txt")
        write_stats(usa_bam, out / "map_to_USA300.samtools_stats.txt")
        usa_fs = (out / "map_to_USA300.flagstat.txt").read_text()
        tot, mp, pct = parse_flagstat_mapped_pct(usa_fs)

        write_flagstat(pan_bam, out / "map_to_panel.flagstat.txt")
        panel_hit_counts(pan_bam, out / "map_to_panel.primary_hits.tsv")

        cov_raw = out / "map_to_contigs.samtools_coverage.tsv"
        samtools_coverage(ctg_bam, cov_raw)
        merge_gc_coverage(args.assembly, cov_raw, out / "contigs_gc_depth_merged.tsv")

        # Summary line
        panel_counts = {}
        with (out / "map_to_panel.primary_hits.tsv").open() as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            for row in rd:
                panel_counts[row["panel_contig"]] = int(row["primary_alignments"])
        staph_n = panel_counts.get("Staphylococcus_aureus_CP000046", 0)
        ecoli_n = panel_counts.get("Escherichia_coli_NC000913", 0)
        panel_sum = staph_n + ecoli_n

        with (out / "diagnose_summary.txt").open("w") as fh:
            fh.write(f"sample\t{args.sample}\n")
            fh.write(f"platform\t{plat}\n")
            fh.write(f"map_to_USA300_total_reads\t{tot}\n")
            fh.write(f"map_to_USA300_mapped_reads\t{mp}\n")
            fh.write(f"map_to_USA300_mapped_pct\t{pct}\n")
            fh.write(f"panel_staph_primary_hits\t{staph_n}\n")
            fh.write(f"panel_ecoli_primary_hits\t{ecoli_n}\n")
            fh.write(f"panel_assigned_sum\t{panel_sum}\n")
            if panel_sum:
                fh.write(
                    f"panel_frac_staph\t{100.0 * staph_n / panel_sum:.2f}\n"
                )
                fh.write(
                    f"panel_frac_ecoli\t{100.0 * ecoli_n / panel_sum:.2f}\n"
                )

        pct_s = f"{pct:.1f}" if pct is not None else "NA"
        print(
            f"{args.sample}: USA300 mapped {pct_s}% | "
            f"panel Staph/Ecoli hits {staph_n}/{ecoli_n} → {out}"
        )
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
