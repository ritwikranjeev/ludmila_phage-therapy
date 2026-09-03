#!/usr/bin/env python3
"""Gene presence/absence by reciprocal best hit with a per-cluster score floor.

WHY THIS REPLACES BOTH EARLIER CALLS
------------------------------------
Panaroo splits a gene across clusters when its alleles diverge, so each arm looks
like an accessory gene absent from the strains carrying the other arm. sraP is
split 73/10 that way; 1,529 same-annotation cluster pairs in this pangenome are
complementary, and 30% of "accessory" clusters turned out to be core.

A plain protein search fixes that but overshoots: searching with a tarM query at
E < 1e-20 and 50% coverage also hits UDP-N-acetylglucosamine--peptide
N-acetylglucosaminyltransferase, a different enzyme in the same fold family. So
Panaroo under-calls and a permissive search over-calls.

This script calls a gene present in a strain only when all three hold:

  1. SCORE FLOOR. The hit scores at least `--floor-fraction` of the weakest score
     the query gets against the cluster's own members. Members are known
     positives, so this calibrates the threshold per family rather than using one
     E-value for a 100 aa protein and a 2,400 aa surface protein alike.
  2. COVERAGE. The alignment covers at least `--min-coverage` of the query.
  3. RECIPROCAL BEST HIT. Of all cluster representatives in the pangenome, the
     one this protein matches best is the cluster being tested. A paralogue fails
     here because it matches its own cluster better.

Output matches recall_presence_absence.py: cluster, query length, 0/1 per strain.

USAGE
    python3 recall_rbh.py --host-range /path/to/Host-Range --out work/recall_rbh.tsv --threads 8
"""
import argparse, csv, os, sys, time
import numpy as np

try:
    from pyhmmer.easel import TextSequence, Alphabet, DigitalSequenceBlock
    from pyhmmer.plan7 import Pipeline
except ImportError:
    sys.exit("pyhmmer is required:  pip install pyhmmer")

ABC = Alphabet.amino()
MIN_STRAINS = 5


def read_faa(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def load(host_range):
    pan = os.path.join(host_range, "data/hosts/panaroo_pangenome/pangenome")
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(pan, "gene_presence_absence.Rtab"))]
    strains = rows[0][1:]
    counts = {r[0]: sum(int(x) for x in r[1:]) for r in rows[1:]}
    n = len(strains)
    testable = [g for g, c in counts.items() if MIN_STRAINS <= c <= n - MIN_STRAINS]

    members = {}
    with open(os.path.join(pan, "gene_presence_absence.csv")) as fh:
        rdr = csv.reader(fh)
        hdr = next(rdr)
        for row in rdr:
            ids = {hdr[i]: row[i].split(";")[0] for i in range(3, len(row)) if row[i]}
            if ids:
                members[row[0]] = ids
    return strains, testable, members


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--host-range", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--floor-fraction", type=float, default=0.8,
                    help="fraction of the weakest within-cluster score a hit must reach")
    ap.add_argument("--min-coverage", type=float, default=0.7)
    args = ap.parse_args()

    strains, testable, members = load(args.host_range)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    with open(os.path.join(os.path.dirname(os.path.abspath(args.out)),
                           "recall_strain_order.txt"), "w") as fh:
        fh.write("\n".join(strains))

    print(f"loading {len(strains)} proteomes", flush=True)
    prot = {s: read_faa(os.path.join(args.host_range, "data/hosts/prokka_annotations",
                                     s, f"{s}.faa")) for s in strains}

    # One representative per cluster: the longest member. Used both as the query
    # and, collectively, as the database for the reciprocal check.
    rep_seq, rep_of = {}, {}
    for cl, ids in members.items():
        best, blen = None, 0
        for s, gid in ids.items():
            sq = prot.get(s, {}).get(gid)
            if sq and len(sq) > blen:
                best, blen = sq, len(sq)
        if best:
            rep_seq[cl] = best
    print(f"{len(rep_seq)} cluster representatives, {len(testable)} clusters to call", flush=True)

    # Score floor per cluster: weakest score the representative gets against its
    # own members.
    floors = {}
    for cl in testable:
        if cl not in rep_seq:
            continue
        mem_seqs = [prot[s][g] for s, g in members[cl].items()
                    if g in prot.get(s, {})]
        if not mem_seqs:
            continue
        q = TextSequence(name=cl.encode(), sequence=rep_seq[cl]).digitize(ABC)
        db = DigitalSequenceBlock(ABC, [TextSequence(name=f"m{i}".encode(), sequence=m).digitize(ABC)
                                        for i, m in enumerate(mem_seqs)])
        hits = Pipeline(ABC, E=1e5).search_seq(q, db)
        scores = [h.score for h in hits]
        floors[cl] = min(scores) * args.floor_fraction if scores else 0.0
    print(f"score floors computed for {len(floors)} clusters", flush=True)

    done = set()
    if os.path.exists(args.out):
        done = {l.split("\t")[0] for l in open(args.out)}
    todo = [c for c in testable if c in floors and c not in done]
    print(f"{len(done)} already done, {len(todo)} to go", flush=True)

    # All representatives as one database for the reciprocal direction.
    rep_names = list(rep_seq)
    rep_db = DigitalSequenceBlock(
        ABC, [TextSequence(name=c.encode(), sequence=rep_seq[c]).digitize(ABC) for c in rep_names])

    t0 = time.time()
    with open(args.out, "a") as out:
        for i, cl in enumerate(todo, 1):
            q = TextSequence(name=cl.encode(), sequence=rep_seq[cl]).digitize(ABC)
            qlen = len(rep_seq[cl])
            calls = []
            for s in strains:
                db = DigitalSequenceBlock(
                    ABC, [TextSequence(name=k.encode(), sequence=v).digitize(ABC)
                          for k, v in prot[s].items()])
                hits = Pipeline(ABC, E=1e-5).search_seq(q, db)
                present = 0
                if len(hits):
                    h = hits[0]
                    cov = sum(d.alignment.target_to - d.alignment.target_from + 1
                              for d in h.domains.included) / qlen
                    if h.score >= floors[cl] and cov >= args.min_coverage:
                        # reciprocal check: does this protein match this cluster best?
                        hid = h.name.decode() if isinstance(h.name, bytes) else str(h.name)
                        back = Pipeline(ABC, E=1e-5).search_seq(
                            TextSequence(name=b"q", sequence=prot[s][hid]).digitize(ABC), rep_db)
                        if len(back):
                            top = back[0].name
                            top = top.decode() if isinstance(top, bytes) else str(top)
                            present = int(top == cl)
                calls.append(present)
            out.write(f"{cl}\t{qlen}\t{''.join(map(str, calls))}\n")
            out.flush()
            if i % 25 == 0:
                rate = (time.time() - t0) / i
                print(f"  {i}/{len(todo)}  {rate:.1f}s per cluster  "
                      f"~{rate*(len(todo)-i)/60:.0f} min left", flush=True)
    print(f"done in {(time.time()-t0)/60:.1f} min -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
