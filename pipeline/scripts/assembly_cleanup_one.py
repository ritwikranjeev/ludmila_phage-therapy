#!/usr/bin/env python3
"""
Per-sample assembly cleanup: drop contigs that are neither well-supported by reads
nor align to a *S. aureus* USA300 chromosome reference (CP000046.1).

Keeps likely plasmids if read depth is comparable to the chromosomal median.
"""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from statistics import median
from typing import Dict, List, Optional, Tuple
from urllib.request import urlopen


def read_fasta(path: Path) -> Dict[str, str]:
    seqs: Dict[str, List[str]] = {}
    cur: Optional[str] = None
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs[cur] = []
            elif cur:
                seqs[cur].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def write_fasta(path: Path, seqs: Dict[str, str], order: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        for cid in order:
            if cid not in seqs:
                continue
            seq = seqs[cid]
            out.write(f">{cid}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i : i + 60] + "\n")


def ensure_reference(ref_path: Path) -> None:
    if ref_path.is_file() and ref_path.stat().st_size > 1000:
        return
    ref_path.parent.mkdir(parents=True, exist_ok=True)
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        "db=nucleotide&id=CP000046.1&rettype=fasta&retmode=text"
    )
    print(f"Downloading reference {ref_path.name} from NCBI …", file=sys.stderr)
    with urlopen(url, timeout=120) as resp:
        data = resp.read().decode()
    ref_path.write_text(data)
    if ref_path.stat().st_size < 1000:
        sys.exit("ERROR: reference download failed or too small")


def run(cmd: List[str], cwd: Optional[Path] = None) -> None:
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write(r.stdout or "")
        sys.stderr.write(r.stderr or "")
        sys.exit(f"ERROR: command failed ({r.returncode}): {' '.join(cmd)}")


def minimap2_contigs_vs_ref(
    ref: Path, contigs: Path, paf_out: Path, threads: int
) -> None:
    run(
        [
            "minimap2",
            "-t",
            str(threads),
            "-c",
            "--secondary=no",
            "-x",
            "asm10",
            str(ref),
            str(contigs),
            "-o",
            str(paf_out),
        ]
    )


def parse_paf_best(paf_path: Path, min_id: float, min_frac: float) -> Dict[str, bool]:
    best: Dict[str, Tuple[float, float]] = {}
    if not paf_path.is_file() or paf_path.stat().st_size == 0:
        return {}
    with paf_path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            qname, qlen, qstart, qend = p[0], int(p[1]), int(p[2]), int(p[3])
            matches, alen = int(p[9]), int(p[10])
            if alen <= 0 or qlen <= 0:
                continue
            ident = matches / alen
            span = (qend - qstart) / qlen
            frac = max(span, alen / qlen)
            if qname not in best or ident > best[qname][0]:
                best[qname] = (ident, frac)
    return {q: (ident >= min_id and frac >= min_frac) for q, (ident, frac) in best.items()}


def map_reads_long(
    work_contigs: Path, reads: Path, bam_out: Path, threads: int, tmpdir: Path
) -> None:
    sam = tmpdir / "long.sam"
    p1 = subprocess.run(
        [
            "minimap2",
            "-ax",
            "map-ont",
            "-t",
            str(threads),
            str(work_contigs),
            str(reads),
        ],
        capture_output=True,
    )
    if p1.returncode != 0:
        sys.stderr.write(p1.stderr.decode())
        sys.exit("ERROR: minimap2 long read mapping failed")
    sam.write_bytes(p1.stdout)
    run(
        [
            "samtools",
            "sort",
            "-@",
            str(threads),
            "-o",
            str(bam_out),
            str(sam),
        ]
    )
    sam.unlink(missing_ok=True)
    run(["samtools", "index", str(bam_out)])


def map_reads_illumina(
    work_contigs: Path, r1: Path, r2: Path, bam_out: Path, threads: int, tmpdir: Path
) -> None:
    run(["bwa", "index", str(work_contigs)], cwd=tmpdir)
    sam = tmpdir / "pe.sam"
    p1 = subprocess.run(
        [
            "bwa",
            "mem",
            "-t",
            str(threads),
            str(work_contigs.name),
            str(r1.resolve()),
            str(r2.resolve()),
        ],
        cwd=tmpdir,
        capture_output=True,
    )
    if p1.returncode != 0:
        sys.stderr.write(p1.stderr.decode())
        sys.exit("ERROR: bwa mem failed")
    sam.write_bytes(p1.stdout)
    run(
        [
            "samtools",
            "sort",
            "-@",
            str(threads),
            "-o",
            str(bam_out),
            str(sam),
        ]
    )
    sam.unlink(missing_ok=True)
    run(["samtools", "index", str(bam_out)])


def samtools_coverage_mean_depth(bam: Path, cov_tsv: Path) -> Dict[str, float]:
    run(
        [
            "samtools",
            "coverage",
            "-o",
            str(cov_tsv),
            str(bam),
        ]
    )
    depths: Dict[str, float] = {}
    with cov_tsv.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            rname, meandepth = parts[0], float(parts[6])
            depths[rname] = meandepth
    return depths


def rm_bwa_index(stem: Path) -> None:
    for suf in (".amb", ".ann", ".bwt", ".pac", ".sa"):
        p = Path(str(stem) + suf)
        p.unlink(missing_ok=True)


def load_reads_from_samples_tsv(samples_tsv: Path, sample: str) -> Tuple[Optional[Path], Optional[Path], Optional[Path]]:
    with samples_tsv.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("sample") != sample:
                continue
            r1, r2, lr = row.get("r1", ""), row.get("r2", ""), row.get("longread", "")
            p1 = Path(r1) if r1 and r1 != "NA" else None
            p2 = Path(r2) if r2 and r2 != "NA" else None
            plr = Path(lr) if lr and lr != "NA" else None
            return p1, p2, plr
    sys.exit(f"ERROR: sample {sample!r} not found in {samples_tsv}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--contigs", type=Path, required=True)
    ap.add_argument("--samples-tsv", type=Path, default=None, help="If set, read r1/r2/longread paths from this TSV")
    ap.add_argument("--r1", type=Path, default=None)
    ap.add_argument("--r2", type=Path, default=None)
    ap.add_argument("--longread", type=Path, default=None)
    ap.add_argument("--ref", type=Path, required=True)
    ap.add_argument("--out-fasta", type=Path, required=True)
    ap.add_argument("--out-removed", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--min-identity", type=float, default=0.88)
    ap.add_argument("--min-frac-aln", type=float, default=0.30)
    ap.add_argument("--cov-rel", type=float, default=0.22)
    ap.add_argument("--min-abs-depth", type=float, default=1.5)
    ap.add_argument("--large-keep-bp", type=int, default=80_000)
    args = ap.parse_args()

    if args.samples_tsv:
        r1, r2, lr = load_reads_from_samples_tsv(args.samples_tsv, args.sample)
        args.r1, args.r2, args.longread = r1, r2, lr

    ensure_reference(args.ref)

    contigs = read_fasta(args.contigs)
    if not contigs:
        sys.exit("ERROR: empty contigs")

    tmpdir = Path(tempfile.mkdtemp(prefix="asmclean_"))
    work_ctg = tmpdir / "contigs.fa"
    try:
        shutil.copy2(args.contigs, work_ctg)

        paf = tmpdir / "ctg_vs_ref.paf"
        minimap2_contigs_vs_ref(args.ref, work_ctg, paf, args.threads)
        ref_ok = parse_paf_best(paf, args.min_identity, args.min_frac_aln)

        bam = tmpdir / "reads_to_ctg.bam"
        lr = args.longread
        if lr and str(lr) != "NA" and lr.is_file():
            map_reads_long(work_ctg, lr, bam, args.threads, tmpdir)
        elif (
            args.r1
            and args.r2
            and args.r1.is_file()
            and args.r2.is_file()
        ):
            map_reads_illumina(work_ctg, args.r1, args.r2, bam, args.threads, tmpdir)
            rm_bwa_index(work_ctg)
        else:
            sys.exit("ERROR: no usable reads (long or paired) for sample")

        cov_tsv = tmpdir / "coverage.tsv"
        depths = samtools_coverage_mean_depth(bam, cov_tsv)

        long_for_median = [
            depths[c]
            for c, s in contigs.items()
            if len(s) >= 5000 and c in depths
        ]
        if not long_for_median:
            long_for_median = list(depths.values())
        med = median(long_for_median) if long_for_median else 0.0
        if med <= 0:
            med = max(depths.values()) if depths else 1.0

        keep_order: List[str] = []
        removed_rows: List[Tuple[str, str, int, float, str]] = []

        for cid, seq in contigs.items():
            ln = len(seq)
            d = depths.get(cid, 0.0)
            hits_ref = ref_ok.get(cid, False)
            # OR not AND: requiring max(min_abs, frac*median) forces short contigs to match
            # the full median depth and wipes fragmented Illumina assemblies (e.g. 8k contigs).
            hi_cov = (d >= args.min_abs_depth) or (
                med > 0 and d >= args.cov_rel * med
            )
            large = ln >= args.large_keep_bp

            if hits_ref:
                keep_order.append(cid)
            elif hi_cov:
                keep_order.append(cid)
            elif large:
                keep_order.append(cid)
            else:
                if d < args.min_abs_depth and not hits_ref:
                    reason = "drop_low_depth"
                else:
                    reason = "drop_depth_below_median_fraction"
                removed_rows.append((args.sample, cid, ln, d, reason))

        write_fasta(args.out_fasta, contigs, keep_order)

        args.out_removed.parent.mkdir(parents=True, exist_ok=True)
        with args.out_removed.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["sample", "contig", "length", "meandepth", "reason"])
            for row in removed_rows:
                w.writerow(row)

        print(
            f"{args.sample}: kept {len(keep_order)}/{len(contigs)} contigs; "
            f"removed {len(removed_rows)} (median depth ~{med:.2f})"
        )
    finally:
        for p in sorted(tmpdir.glob("*"), reverse=True):
            try:
                if p.is_file():
                    p.unlink()
            except OSError:
                pass
        try:
            tmpdir.rmdir()
        except OSError:
            pass


if __name__ == "__main__":
    main()
