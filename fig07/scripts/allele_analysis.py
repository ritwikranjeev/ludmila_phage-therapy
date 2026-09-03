#!/usr/bin/env python3
"""Allele-level association for the WTA/LTA receptor pathway.

WHY ALLELES AND NOT GENE CONTENT
--------------------------------
Once genes are called properly (recall_genes.py), the receptor pathway turns out
to be almost invariant in gene content: tarS, tarL, tagH, fmtA, ltaS and agrB are
present in 81-83 of 83 strains. Only tarM is genuinely variable (33/83). So a
presence/absence scan cannot address the paper's hypothesis, because there is
nothing to vary. If wall teichoic acid mediates adsorption, the question is which
ALLELE a strain carries.

WHAT THIS TESTS, per gene in the pathway:
  * every variable amino-acid position (minor allele count >= --min-count)
  * truncation: protein shorter than --trunc-fraction of the gene's median length,
    which is how a frameshift or nonsense variant shows up in a proteome
  * divergence: identity to the consensus allele, as a quantitative feature

All are tested with the same mixed model as the genome-wide scan (core-genome
kinship, genome size and CDS count as covariates), and the family-wise threshold
comes from permuting the phenotype across the whole pre-specified family — so the
multiplicity of testing every position in every gene is paid for honestly.

USAGE
    python3 allele_analysis.py --host-range /path/to/Host-Range --work work \\
        [--genes tarM,tarS,tarP,tarK,tarL,tagH,tagG,tagX,tagU,fmtA,ltaS,ltaA,ugtP,mnaA] \\
        [--permutations 500]
"""
import argparse, csv, os, sys, collections
import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from association import (core_genome, assembly_metrics, scan, bh, lambda_gc,
                         permutation_threshold)

try:
    from pyhmmer.easel import TextSequence, Alphabet, DigitalSequenceBlock, MSA, TextMSA
    from pyhmmer.plan7 import Builder, Background
    from pyhmmer import hmmer
except ImportError:
    sys.exit("pyhmmer is required:  pip install pyhmmer")

ABC = Alphabet.amino()
DEFAULT_GENES = ["tarM", "tarS", "tarP", "tarK", "tarL", "tarA", "tarB", "tarD",
                 "tarF", "tagH", "tagG", "tagX", "tagU", "fmtA", "ltaS", "ltaA",
                 "ugtP", "mnaA", "pgcA", "gtaB", "sraP", "clfA", "clfB", "sdrE", "spa"]


def as_str(x):
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


def gene_sequences(host_range, strains, gene_names):
    """Per-strain protein for each requested gene, via Panaroo cluster membership."""
    pan = os.path.join(host_range, "data/hosts/panaroo_pangenome/pangenome")
    members = {}
    with open(os.path.join(pan, "gene_presence_absence.csv")) as fh:
        rdr = csv.reader(fh)
        hdr = next(rdr)
        for row in rdr:
            base = {t.split("_")[0] for t in row[0].split("~~~")}
            hit = base & set(gene_names)
            if not hit:
                continue
            g = sorted(hit)[0]
            for i in range(3, len(row)):
                if row[i]:
                    members.setdefault(g, {}).setdefault(hdr[i], row[i].split(";")[0])
    out = {g: {} for g in members}
    for s in strains:
        prot = read_faa(os.path.join(host_range, "data/hosts/prokka_annotations", s, f"{s}.faa"))
        for g, ids in members.items():
            gid = ids.get(s)
            if gid and gid in prot:
                out[g][s] = prot[gid]
    return out


def align(seqs_by_strain):
    """Align a gene's alleles by building an HMM from the longest allele."""
    strains = list(seqs_by_strain)
    longest = max(strains, key=lambda s: len(seqs_by_strain[s]))
    seed = TextMSA(name=b"seed",
                   sequences=[TextSequence(name=longest.encode(),
                                           sequence=seqs_by_strain[longest])]).digitize(ABC)
    hmm, _, _ = Builder(ABC).build_msa(seed, Background(ABC))
    digital = [TextSequence(name=s.encode(), sequence=seqs_by_strain[s]).digitize(ABC)
               for s in strains]
    msa = hmmer.hmmalign(hmm, digital)
    text = msa if isinstance(msa, TextMSA) else msa.textize()
    aln = {as_str(text.names[i]): str(text.alignment[i]).upper()
           for i in range(len(text.names))}
    return aln


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--host-range", required=True)
    ap.add_argument("--work", default="work")
    ap.add_argument("--genes", default=",".join(DEFAULT_GENES))
    ap.add_argument("--min-count", type=int, default=5)
    ap.add_argument("--trunc-fraction", type=float, default=0.9)
    ap.add_argument("--permutations", type=int, default=500)
    args = ap.parse_args()
    W = args.work
    genes = [g.strip() for g in args.genes.split(",") if g.strip()]

    pheno = list(csv.DictReader(
        open(os.path.join(args.host_range, "config/meta/host-range-results.tsv")),
        delimiter="\t"))
    strains = [r["sample"] for r in pheno]
    y = np.array([int(r["infection"]) for r in pheno], float)
    K, D = core_genome(args.host_range, strains,
                       os.path.join(W, "kinship.npy"), os.path.join(W, "snp_distance.npy"))
    asm = assembly_metrics(args.host_range, strains, os.path.join(W, "assembly_metrics.tsv"))
    z = lambda v: (v - v.mean()) / v.std()
    cov = np.column_stack([z(asm[:, 1]), z(asm[:, 3])])

    print(f"extracting alleles for {len(genes)} genes", flush=True)
    seqs = gene_sequences(args.host_range, strains, genes)
    print("gene            strains  length range", flush=True)
    for g in sorted(seqs):
        L = [len(v) for v in seqs[g].values()]
        print(f"  {g:12s} {len(seqs[g]):4d}     {min(L)}-{max(L)} aa", flush=True)

    features, labels = [], []
    for g in sorted(seqs):
        by_strain = seqs[g]
        if len(by_strain) < 20:
            continue
        lengths = np.array([len(by_strain.get(s, "")) for s in strains], float)
        med = np.median([l for l in lengths if l > 0])

        # truncation, as a proxy for frameshift or nonsense
        trunc = ((lengths > 0) & (lengths < args.trunc_fraction * med)).astype(float)
        if args.min_count <= trunc.sum() <= len(strains) - args.min_count:
            features.append(trunc)
            labels.append(f"{g}:truncated")

        aln = align(by_strain)
        cols = len(next(iter(aln.values())))
        mat = np.array([list(aln.get(s, "-" * cols)) for s in strains])
        # per-position amino-acid variants
        for j in range(cols):
            col = mat[:, j]
            for aa in set(col) - {"-", "."}:
                v = (col == aa).astype(float)
                if args.min_count <= v.sum() <= len(strains) - args.min_count:
                    features.append(v)
                    labels.append(f"{g}:pos{j+1}:{aa}")
        # divergence from the consensus allele
        consensus = np.array([collections.Counter(mat[:, j]).most_common(1)[0][0]
                              for j in range(cols)])
        ident = np.array([np.mean(mat[i] == consensus) for i in range(len(strains))])
        if ident.std() > 0:
            features.append(z(ident))
            labels.append(f"{g}:identity_to_consensus")

    M = np.array(features)
    print(f"\n{len(labels)} allele features across {len(seqs)} genes", flush=True)
    delta, res = scan(y, K, cov, M, labels, min_count=0)
    ps = np.array([r[3] for r in res])
    q = bh(ps)
    print(f"tested={len(res)}  lambda={lambda_gc(ps):.2f}", flush=True)

    mins = permutation_threshold(y, K, cov, M, labels, args.permutations)
    thr = mins[int(0.05 * len(mins))]
    print(f"{args.permutations} permutations -> family-wise 5% threshold p = {thr:.2e}")

    with open(os.path.join(W, "allele_association.tsv"), "w") as fh:
        fh.write("feature\tn_present\tbeta\tp\tq_bh\tfwer_p\n")
        for (lab, c, b, p), qq in zip(res, q):
            fh.write(f"{lab}\t{c}\t{b:.4f}\t{p:.3e}\t{qq:.4f}\t{(mins <= p).mean():.4f}\n")
    o = np.argsort(ps)
    print(f"\n{'feature':44s}{'n+':>5s}{'beta':>8s}{'p':>10s}{'FWER':>8s}")
    for i in o[:12]:
        lab, c, b, p = res[i]
        print(f"{lab[:44]:44s}{c:5d}{b:8.3f}{p:10.2e}{(mins <= p).mean():8.3f}")
    print(f"\nsignificant at family-wise 5%: {(ps <= thr).sum()}")


if __name__ == "__main__":
    main()
