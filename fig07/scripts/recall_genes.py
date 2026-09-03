#!/usr/bin/env python3
"""Gene presence/absence, calling genes rather than Panaroo clusters.

THE PROBLEM THIS SOLVES
-----------------------
Three ways of calling the same thing disagree, and each fails differently:

  Panaroo clusters      splits a gene when its alleles diverge, so each arm looks
                        accessory. sraP is split 73/10; rplI is split 24/59 and
                        the two arms have OPPOSITE associations with phenotype.
  Permissive homology   fixes the splits but merges paralogues: a tarM query also
                        hits UDP-N-acetylglucosamine--peptide transferase, giving
                        the nonsense answer that tarM is in 83/83 strains.
  Reciprocal best hit   fixes paralogues but re-breaks splits: a strain carrying
                        arm B of sraP fails the reciprocal test for arm A, and
                        sraP collapses to 7/83.

So the calling has to happen at the level of the GENE, not the cluster. This
script therefore runs in two stages:

  STAGE 1  Merge cluster representatives into gene families, by two rules:
             (a) mutual best hit, and
             (b) sequence-related (best hit, >=70% coverage) AND never co-occurring
                 in the same strain, which is the signature of one gene split in
                 two rather than two genes.
           Rule (b) is what recovers rplI, whose two arms are 24 and 59 strains
           with 2 co-occurrences; rule (a) alone missed it. tarM and tarS stay
           separate under both rules, because they co-occur freely.
  STAGE 2  For each strain, assign every protein to the family it matches best,
           then call a family present if some protein is assigned to it and
           clears that family's score floor and coverage.

The score floor is calibrated per family from its own members — the weakest
score the representative achieves against known members — so a 100 aa protein
and a 2,400 aa surface protein are not judged by the same E-value.

USAGE
    python3 recall_genes.py --host-range /path/to/Host-Range --out work/genes.tsv --threads 8

Outputs `--out` (family, representative cluster, members, query length, 0/1 per
strain) and, alongside it, `gene_families.tsv` recording which clusters were
merged into which family.
"""
import argparse, csv, os, sys, time, collections
import numpy as np

try:
    from pyhmmer.easel import TextSequence, Alphabet, DigitalSequenceBlock
    from pyhmmer.plan7 import Pipeline
    from pyhmmer import hmmer
except ImportError:
    sys.exit("pyhmmer is required:  pip install pyhmmer")

ABC = Alphabet.amino()


def as_str(x):
    """pyhmmer returns names as bytes or str depending on version/path."""
    return x.decode() if isinstance(x, (bytes, bytearray)) else str(x)

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


def digital(seqs):
    return [TextSequence(name=k.encode(), sequence=v).digitize(ABC) for k, v in seqs.items()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--host-range", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--floor-fraction", type=float, default=0.8)
    ap.add_argument("--min-coverage", type=float, default=0.7)
    ap.add_argument("--max-families", type=int, default=0,
                    help="debug: only call this many families (0 = all)")
    ap.add_argument("--query-chunk", type=int, default=500,
                    help="family queries per search batch; lower it if memory is tight")
    args = ap.parse_args()

    pan = os.path.join(args.host_range, "data/hosts/panaroo_pangenome/pangenome")
    outdir = os.path.dirname(os.path.abspath(args.out)) or "."
    os.makedirs(outdir, exist_ok=True)

    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(pan, "gene_presence_absence.Rtab"))]
    strains = rows[0][1:]
    with open(os.path.join(outdir, "recall_strain_order.txt"), "w") as fh:
        fh.write("\n".join(strains))

    members = {}
    with open(os.path.join(pan, "gene_presence_absence.csv")) as fh:
        rdr = csv.reader(fh)
        hdr = next(rdr)
        for row in rdr:
            ids = {hdr[i]: row[i].split(";")[0] for i in range(3, len(row)) if row[i]}
            if ids:
                members[row[0]] = ids

    # One streaming pass over the proteomes: keep only the longest member of each
    # cluster (the representative) and up to MAX_FLOOR_MEMBERS others for
    # calibrating that family's score floor. Holding all 83 proteomes at once
    # costs a few GB and gets the process OOM-killed on a small machine.
    MAX_FLOOR_MEMBERS = 5
    print(f"scanning {len(strains)} proteomes", flush=True)
    rep, rep_len, floor_members = {}, {}, collections.defaultdict(list)
    for s in strains:
        seqs = read_faa(os.path.join(args.host_range, "data/hosts/prokka_annotations",
                                     s, f"{s}.faa"))
        for cl, ids in members.items():
            gid = ids.get(s)
            if not gid:
                continue
            sq = seqs.get(gid)
            if not sq:
                continue
            if len(sq) > rep_len.get(cl, 0):
                rep[cl], rep_len[cl] = sq, len(sq)
            if len(floor_members[cl]) < MAX_FLOOR_MEMBERS:
                floor_members[cl].append(sq)
        del seqs
    rep_names = list(rep)
    print(f"{len(rep_names)} cluster representatives", flush=True)

    # ---- stage 1: merge clusters into gene families by mutual best hit -------
    t0 = time.time()
    rep_db = [TextSequence(name=c.encode(), sequence=rep[c]).digitize(ABC) for c in rep_names]
    queries = [TextSequence(name=c.encode(), sequence=rep[c]).digitize(ABC) for c in rep_names]
    best_of = {}
    for hits in hmmer.phmmer(queries, rep_db, cpus=args.threads, E=1e-5):
        qname = as_str(hits.query.name)
        for h in hits:
            hname = as_str(h.name)
            if hname == qname:
                continue
            cov = sum(d.alignment.target_to - d.alignment.target_from + 1
                      for d in h.domains.included) / len(rep[qname])
            if cov >= args.min_coverage:
                best_of[qname] = hname
            break
    parent = {c: c for c in rep_names}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
            return 1
        return 0

    rule1 = sum(union(a, b) for a, b in best_of.items() if best_of.get(b) == a)

    # rule (b): related, but never present in the same strain -> one gene split
    pres = {r[0]: np.array([int(x) for x in r[1:]], bool) for r in rows[1:]}
    rule2 = 0
    for a, b in best_of.items():
        if a not in pres or b not in pres:
            continue
        inter = int((pres[a] & pres[b]).sum())
        uni = int((pres[a] | pres[b]).sum())
        if inter <= 2 and uni >= 0.5 * len(strains):
            rule2 += union(a, b)

    fam = collections.defaultdict(list)
    for c in rep_names:
        fam[find(c)].append(c)
    print(f"stage 1: {rule1} mutual-best-hit merges, {rule2} split-arm merges "
          f"-> {len(fam)} gene families ({time.time()-t0:.0f}s)", flush=True)

    if args.max_families:
        probe = [f for f in fam
                 if any(c.startswith(("sraP", "tarM", "tarS", "rplI", "agrB")) for c in fam[f])]
        keep = list(dict.fromkeys(probe + list(fam)))[:args.max_families]
        fam = {f: fam[f] for f in keep}
        keep_clusters = {c for cls in fam.values() for c in cls}
        floor_members = {c: v for c, v in floor_members.items() if c in keep_clusters}
        print(f"debug: restricted to {len(fam)} families", flush=True)

    with open(os.path.join(outdir, "gene_families.tsv"), "w") as fh:
        fh.write("family\tn_clusters\tclusters\n")
        for f, cls in sorted(fam.items(), key=lambda kv: -len(kv[1])):
            fh.write(f"{f}\t{len(cls)}\t{','.join(sorted(cls))}\n")

    # family representative = longest member across its clusters
    fam_rep = {}
    for f, cls in fam.items():
        best, blen = None, 0
        for cl in cls:
            if len(rep[cl]) > blen:
                best, blen = rep[cl], len(rep[cl])
        fam_rep[f] = best

    # Score floor per family, computed in ONE batched search rather than 4,372:
    # every sampled member is put in a single database, named "<cluster>|<i>", and
    # each family's floor is the weakest score its representative achieves against
    # members of its own clusters.
    print("building family representatives and score floors", flush=True)
    cluster_family = {cl: f for f, cls in fam.items() for cl in cls}
    member_db, member_owner = [], {}
    for cl, seqs in floor_members.items():
        for i, m in enumerate(seqs):
            nm = f"{cl}|{i}"
            member_owner[nm] = cluster_family.get(cl)
            member_db.append(TextSequence(name=nm.encode(), sequence=m).digitize(ABC))
    fam_names_all = list(fam)
    fam_floor = {f: None for f in fam_names_all}
    t1 = time.time()
    for i in range(0, len(fam_names_all), args.query_chunk):
        batch = [TextSequence(name=f.encode(), sequence=fam_rep[f]).digitize(ABC)
                 for f in fam_names_all[i:i + args.query_chunk]]
        for hits in hmmer.phmmer(batch, member_db, cpus=args.threads, E=1e5):
            f = as_str(hits.query.name)
            own = [h.score for h in hits if member_owner.get(as_str(h.name)) == f]
            if own:
                fam_floor[f] = min(own) * args.floor_fraction
    for f in fam_names_all:
        if fam_floor[f] is None:
            fam_floor[f] = 0.0
    print(f"score floors: {len(fam_floor)} families in {time.time()-t1:.0f}s", flush=True)
    del floor_members, rep, member_db


    # ---- stage 2: assign each strain's proteins to families -----------------
    fam_names = list(fam)
    fam_queries = [TextSequence(name=f.encode(), sequence=fam_rep[f]).digitize(ABC)
                   for f in fam_names]
    calls = {f: [] for f in fam_names}
    partial = args.out + ".partial"
    done_strains = []
    if os.path.exists(partial):
        for line in open(partial):
            st, bits = line.rstrip("\n").split("\t")
            done_strains.append(st)
            for f, b in zip(fam_names, bits):
                calls[f].append(int(b))
        print(f"resuming: {len(done_strains)} strains already called", flush=True)
    t0 = time.time()
    for si, s in enumerate(strains, 1):
        if s in done_strains:
            continue
        strain_prot = read_faa(os.path.join(args.host_range, "data/hosts/prokka_annotations",
                                            s, f"{s}.faa"))
        sdb = digital(strain_prot)
        # best family per protein, and best protein per family, from one pass
        best_fam_for_prot, hit_for_fam = {}, {}
        chunk = args.query_chunk
        batches = (fam_queries[i:i+chunk] for i in range(0, len(fam_queries), chunk))
        for batch in batches:
          for hits in hmmer.phmmer(batch, sdb, cpus=args.threads, E=1e-5):
            f = as_str(hits.query.name)
            for h in hits:
                pid = as_str(h.name)
                cov = sum(d.alignment.target_to - d.alignment.target_from + 1
                          for d in h.domains.included) / len(fam_rep[f])
                if h.score < fam_floor[f] or cov < args.min_coverage:
                    continue
                if pid not in best_fam_for_prot or h.score > best_fam_for_prot[pid][1]:
                    best_fam_for_prot[pid] = (f, h.score)
                if f not in hit_for_fam or h.score > hit_for_fam[f][1]:
                    hit_for_fam[f] = (pid, h.score)
                break
        bits = []
        for f in fam_names:
            pid = hit_for_fam.get(f, (None, 0))[0]
            v = int(pid is not None and best_fam_for_prot.get(pid, ("", 0))[0] == f)
            calls[f].append(v)
            bits.append(str(v))
        with open(partial, "a") as pf:
            pf.write(s + "\t" + "".join(bits) + "\n")
        del strain_prot, sdb
        print(f"  strain {si}/{len(strains)}  {s}  "
              f"{sum(calls[f][-1] for f in fam_names)} families present", flush=True)

    with open(args.out, "w") as out:
        out.write("family\tn_clusters\tquery_len\tcalls\n")
        for f in fam_names:
            out.write(f"{f}\t{len(fam[f])}\t{len(fam_rep[f])}\t{''.join(map(str, calls[f]))}\n")
    n = len(strains)
    counts = np.array([sum(calls[f]) for f in fam_names])
    print(f"\ndone in {(time.time()-t0)/60:.1f} min -> {args.out}")
    print(f"  core (>=95%): {(counts >= 0.95*n).sum()} | "
          f"testable (5..{n-5}): {((counts >= 5) & (counts <= n-5)).sum()} | "
          f"rare (<5): {(counts < 5).sum()}")


if __name__ == "__main__":
    main()
