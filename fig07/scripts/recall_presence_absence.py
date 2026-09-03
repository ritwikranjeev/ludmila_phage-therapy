#!/usr/bin/env python3
"""Re-call gene presence/absence by protein search, replacing Panaroo's clusters.

WHY THIS EXISTS
---------------
Panaroo assigns a gene to a cluster by sequence similarity. When alleles of one
gene diverge past its threshold, the gene is split across two clusters, and each
arm then looks like an accessory gene that is "absent" from the strains carrying
the other arm. An audit of 25 randomly chosen testable clusters found that 68%
are present in more strains than Panaroo reports, and 24% are actually core
genes (>=95% of strains) split by divergence. sraP is one of them: two clusters
of 73 and 10 strains, which a protein search resolves to 76 of 83.

Associations computed on the raw matrix therefore test allelic divergence while
appearing to test gene content. This script re-calls every testable cluster by
searching the actual proteomes.

METHOD
------
For each cluster: take its longest member protein as the query, search all 83
Prokka proteomes with phmmer, and call the gene present when E < 1e-20 and the
aligned region covers more than half the query.

Output: one line per cluster,   cluster <TAB> query_length_aa <TAB> 0/1 string
with one character per strain, in the order given by --strain-order (written on
first run). Appends as it goes, so re-running resumes where it stopped.

USAGE
-----
    pip install pyhmmer
    python3 recall_presence_absence.py \
        --host-range /path/to/OCPT-Ludmila_Phage-Therapy/Host-Range \
        --out work/recall_all.tsv \
        --threads 8

On a laptop with 8 cores this takes roughly 10-15 minutes; on a cluster node,
a couple of minutes. It is resumable and safe to interrupt.
"""
import argparse, csv, os, sys, time
from multiprocessing import Pool

try:
    from pyhmmer.easel import TextSequence, Alphabet, DigitalSequenceBlock
    from pyhmmer.plan7 import Pipeline
except ImportError:
    sys.exit("pyhmmer is required:  pip install pyhmmer")

ABC = Alphabet.amino()
E_VALUE = 1e-20      # hit must be at least this significant
MIN_COVERAGE = 0.5   # aligned fraction of the query
MIN_STRAINS = 5      # only clusters testable in an association analysis


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def load_inputs(host_range):
    pan = os.path.join(host_range, "data/hosts/panaroo_pangenome/pangenome")
    rtab = os.path.join(pan, "gene_presence_absence.Rtab")
    rows = [l.rstrip("\n").split("\t") for l in open(rtab)]
    strains = rows[0][1:]
    counts = {r[0]: sum(int(x) for x in r[1:]) for r in rows[1:]}
    n = len(strains)
    testable = [g for g, c in counts.items() if MIN_STRAINS <= c <= n - MIN_STRAINS]

    members = {}
    with open(os.path.join(pan, "gene_presence_absence.csv")) as fh:
        rdr = csv.reader(fh)
        hdr = next(rdr)
        csv_strains = hdr[3:]
        want = set(testable)
        for row in rdr:
            if row[0] in want:
                members[row[0]] = {hdr[i]: row[i].split(";")[0]
                                   for i in range(3, len(row)) if row[i]}
    return csv_strains, testable, members


PROT = {}


def init_worker(host_range, strains):
    for s in strains:
        PROT[s] = read_fasta(os.path.join(
            host_range, "data/hosts/prokka_annotations", s, f"{s}.faa"))


def search(args):
    cluster, member_ids, strains = args
    best_seq, best_len = "", 0
    for s, gid in member_ids.items():
        seq = PROT.get(s, {}).get(gid)
        if seq and len(seq) > best_len:
            best_seq, best_len = seq, len(seq)
    if not best_seq:
        return cluster, None
    query = TextSequence(name=cluster.encode(), sequence=best_seq).digitize(ABC)
    calls = []
    for s in strains:
        db = DigitalSequenceBlock(
            ABC, [TextSequence(name=k.encode(), sequence=v).digitize(ABC)
                  for k, v in PROT[s].items()])
        hits = Pipeline(ABC, E=1e-10).search_seq(query, db)
        present = 0
        if len(hits):
            h = hits[0]
            covered = sum(d.alignment.target_to - d.alignment.target_from + 1
                          for d in h.domains.included)
            if h.evalue < E_VALUE and covered / best_len > MIN_COVERAGE:
                present = 1
        calls.append(present)
    return cluster, (best_len, calls)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--host-range", required=True, help="Host-Range directory")
    ap.add_argument("--out", required=True, help="output TSV (appended, resumable)")
    ap.add_argument("--strain-order", default=None,
                    help="where to write the strain order (default: alongside --out)")
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 4)
    args = ap.parse_args()

    strains, testable, members = load_inputs(args.host_range)
    order_path = args.strain_order or os.path.join(
        os.path.dirname(os.path.abspath(args.out)), "recall_strain_order.txt")
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(order_path, "w") as fh:
        fh.write("\n".join(strains))

    done = set()
    if os.path.exists(args.out):
        done = {l.split("\t")[0] for l in open(args.out)}
    todo = [g for g in testable if g not in done]
    print(f"{len(testable)} testable clusters | {len(done)} done | {len(todo)} to go "
          f"| {len(strains)} strains | {args.threads} threads", flush=True)
    if not todo:
        print("nothing to do")
        return

    jobs = [(c, members[c], strains) for c in todo if c in members]
    t0 = time.time()
    with open(args.out, "a") as out, Pool(
            args.threads, initializer=init_worker,
            initargs=(args.host_range, strains)) as pool:
        for i, (cluster, res) in enumerate(
                pool.imap_unordered(search, jobs, chunksize=2), 1):
            if res is None:
                continue
            length, calls = res
            out.write(f"{cluster}\t{length}\t{''.join(map(str, calls))}\n")
            out.flush()
            if i % 50 == 0:
                rate = (time.time() - t0) / i
                print(f"  {i}/{len(jobs)}  {rate:.1f}s per cluster  "
                      f"~{rate * (len(jobs) - i) / 60:.0f} min left", flush=True)
    print(f"done in {(time.time() - t0) / 60:.1f} min -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
