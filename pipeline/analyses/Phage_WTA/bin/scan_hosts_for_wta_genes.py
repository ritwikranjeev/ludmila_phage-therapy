#!/usr/bin/env python3
"""
hmmsearch host Prokka proteomes against the WTA HMM database; write per-hit summaries.

Usage:
    python bin/scan_hosts_for_wta_genes.py              # all hosts -> raw + by_host
    python bin/scan_hosts_for_wta_genes.py --host-id ID # one host (e.g. Slurm array task)
    python bin/scan_hosts_for_wta_genes.py --collect-raw  # merge by_host/* -> raw_wta_hits.tsv

Protein paths: data/host_prokka/<id>/<id>.faa matches Prokka --prefix <id> in --outdir (default
PREFIX.faa). Glob *.faa if the prefix differed.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Use a single pressed combined DB (see build_wta_hmmdb.py). Set False to run
# hmmsearch once per data/wta_hmms/<gene>_ref.hmm file instead (slower, many invocations).
# ---------------------------------------------------------------------------
USE_COMBINED_HMM: bool = True

# HMM hit reporting (full-sequence E-value). Mirror in build_wta_feature_matrix.WTA_EVALUE_MAX.
WTA_EVALUE_MAX: float = 1e-5

def project_root() -> Path:
    return Path(__file__).resolve().parent.parent


def read_host_tsv(tsv_path: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    for line in tsv_path.read_text(encoding="utf-8").splitlines()[1:]:
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t", 1)
        if len(parts) != 2:
            continue
        rows.append((parts[0].strip(), parts[1].strip()))
    return rows


def parse_hmmsearch_domtblout(dom_path: Path) -> dict[str, dict[str, float | int | str]]:
    """
    Parse hmmsearch --domtblout (HMMER 3.4). One best full-sequence hit per query (HMM) name.
    Returns gene -> {evalue, bitscore, hit_length, best_target}.
    """
    best: dict[str, dict[str, float | int | str]] = {}
    text = dom_path.read_text(encoding="utf-8", errors="replace")
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 19:
            continue
        try:
            target = parts[0]
            query = parts[3]
            full_eval = float(parts[6])
            bitscore = float(parts[7])
            alifrom = int(parts[17])
            alito = int(parts[18])
        except (ValueError, IndexError):
            continue
        hit_length = alito - alifrom + 1
        cur = best.get(query)
        if cur is None or full_eval < float(cur["evalue"]):
            best[query] = {
                "evalue": full_eval,
                "bitscore": bitscore,
                "hit_length": hit_length,
                "best_target": target,
            }
    return best


def run_hmmsearch_combined(
    hmm_path: Path,
    protein_fa: Path,
    domtblout: Path,
    cpu: int,
) -> None:
    domtblout.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "hmmsearch",
        "--noali",
        "--cpu",
        str(cpu),
        "--domtblout",
        str(domtblout),
        "-o",
        "/dev/null",
        str(hmm_path),
        str(protein_fa),
    ]
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def run_hmmsearch_individual(
    hmm_dir: Path,
    protein_fa: Path,
    dom_dir: Path,
    cpu: int,
) -> dict[str, dict[str, float | int | str]]:
    """Search each *.hmm except combined; merge best hits."""
    dom_dir.mkdir(parents=True, exist_ok=True)
    merged: dict[str, dict[str, float | int | str]] = {}
    for hmm_file in sorted(hmm_dir.glob("*.hmm")):
        if hmm_file.name.startswith("wta_combined"):
            continue
        domtblout = dom_dir / f"{hmm_file.stem}.domtblout"
        run_hmmsearch_combined(hmm_file, protein_fa, domtblout, cpu)
        for gene, stats in parse_hmmsearch_domtblout(domtblout).items():
            cur = merged.get(gene)
            if cur is None or float(stats["evalue"]) < float(cur["evalue"]):
                merged[gene] = stats
    return merged


def write_hits_tsv(
    path: Path,
    host_id: str,
    best: dict[str, dict[str, float | int | str]],
    *,
    include_header: bool,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    if include_header:
        lines.append("host_id\tgene\tbest_evalue\tbitscore\thit_length\n")
    for gene in sorted(best.keys()):
        st = best[gene]
        ev = st["evalue"]
        evs = f"{ev:.2e}" if isinstance(ev, float) else str(ev)
        lines.append(
            f"{host_id}\t{gene}\t{evs}\t{float(st['bitscore']):.1f}\t{int(st['hit_length'])}\n"
        )
    mode = "w" if include_header else "a"
    with path.open(mode, encoding="utf-8") as fh:
        fh.writelines(lines)


def process_one_host(
    root: Path,
    host_id: str,
    cpu: int,
    *,
    write_by_host: bool,
) -> list[tuple[str, str, str, float, int]]:
    # Prokka output naming may differ; adjust path or glob if needed.
    protein_fa = root / "data" / "host_prokka" / host_id / f"{host_id}.faa"
    if not protein_fa.is_file():
        alts = list((root / "data" / "host_prokka" / host_id).glob("*.faa"))
        if len(alts) == 1:
            protein_fa = alts[0]
        else:
            print(
                f"ERROR: Missing host proteins for {host_id}: expected {protein_fa}",
                file=sys.stderr,
            )
            sys.exit(1)

    hmm_dir = root / "data" / "wta_hmms"
    by_host_dir = root / "data" / "wta_features" / "by_host"
    combined = hmm_dir / "wta_combined.hmm"

    if USE_COMBINED_HMM:
        if not combined.is_file():
            print(f"ERROR: Combined HMM not found: {combined}", file=sys.stderr)
            sys.exit(1)
        domtblout = by_host_dir / f"{host_id}.domtblout"
        run_hmmsearch_combined(combined, protein_fa, domtblout, cpu)
        best = parse_hmmsearch_domtblout(domtblout)
    else:
        dom_sub = by_host_dir / host_id / "per_gene"
        best = run_hmmsearch_individual(hmm_dir, protein_fa, dom_sub, cpu)
        # also write a merged domtblout path for debugging — optional; skip to keep simple

    rows: list[tuple[str, str, str, float, int]] = []
    for gene, st in sorted(best.items()):
        rows.append(
            (
                host_id,
                gene,
                f"{float(st['evalue']):.2e}",
                float(st["bitscore"]),
                int(st["hit_length"]),
            )
        )

    if write_by_host:
        write_hits_tsv(
            by_host_dir / f"{host_id}_wta_hits.tsv",
            host_id,
            best,
            include_header=True,
        )
    return rows


def collect_raw(root: Path) -> None:
    by_host_dir = root / "data" / "wta_features" / "by_host"
    out_path = root / "data" / "wta_features" / "raw_wta_hits.tsv"
    parts = sorted(by_host_dir.glob("*_wta_hits.tsv"))
    if not parts:
        print(f"ERROR: No files matching {by_host_dir}/*_wta_hits.tsv", file=sys.stderr)
        sys.exit(1)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    first = True
    with out_path.open("w", encoding="utf-8") as out:
        for p in parts:
            lines = p.read_text(encoding="utf-8").splitlines()
            if not lines:
                continue
            if first:
                out.write(lines[0] + "\n")
                first = False
                body = lines[1:]
            else:
                body = lines[1:] if lines[0].startswith("host_id") else lines
            for line in body:
                if line.strip():
                    out.write(line + "\n")
    print(f"Wrote {out_path} ({len(parts)} inputs).")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="WTA hmmsearch over host Prokka proteomes.",
        epilog=(
            f"Best hits are written with raw E-values; the feature matrix applies "
            f"WTA_EVALUE_MAX={WTA_EVALUE_MAX} (keep identical to build_wta_feature_matrix.py)."
        ),
    )
    parser.add_argument("--host-id", help="Process a single host ID (Slurm array tasks).")
    parser.add_argument(
        "--collect-raw",
        action="store_true",
        help="Merge data/wta_features/by_host/*_wta_hits.tsv into raw_wta_hits.tsv.",
    )
    args = parser.parse_args()
    root = project_root()

    if args.collect_raw:
        collect_raw(root)
        return 0

    cpu = int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))

    tsv_path = root / "config" / "host_genomes.tsv"
    if not tsv_path.is_file():
        print(f"ERROR: Missing {tsv_path}", file=sys.stderr)
        return 1

    hosts = read_host_tsv(tsv_path)
    if args.host_id:
        hid = args.host_id.strip()
        process_one_host(root, hid, cpu, write_by_host=True)
        return 0

    all_rows: list[tuple[str, str, str, float, int]] = []
    for hid, _fasta in hosts:
        all_rows.extend(process_one_host(root, hid, cpu, write_by_host=True))

    raw_path = root / "data" / "wta_features" / "raw_wta_hits.tsv"
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    with raw_path.open("w", encoding="utf-8") as fh:
        fh.write("host_id\tgene\tbest_evalue\tbitscore\thit_length\n")
        for host_id, gene, evs, bsc, hlen in all_rows:
            fh.write(f"{host_id}\t{gene}\t{evs}\t{bsc:.1f}\t{hlen}\n")
    print(f"Wrote {raw_path} ({len(all_rows)} hit row(s)).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
