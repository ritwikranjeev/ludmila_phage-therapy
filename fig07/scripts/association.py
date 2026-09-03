#!/usr/bin/env python3
"""Association between host gene content and OCPT-Ludmila susceptibility.

Consumes the re-called presence/absence matrix from recall_presence_absence.py
(never Panaroo's raw matrix — see that script for why) and runs:

  1. assembly metrics per strain            -> work/assembly_metrics.tsv
  2. kinship + SNP distances from the core alignment  -> work/kinship.npy, work/snp_distance.npy
  3. genome-wide mixed-model association    -> work/gwas.tsv
  4. pre-specified candidate family         -> work/candidate_family.tsv
  5. within-clone contrasts                 -> work/within_clone.tsv

MODEL
-----
y = mu + covariates + beta * gene + u,   u ~ N(0, sigma_g^2 K),  e ~ N(0, sigma_e^2 I)

K is core-genome relatedness. The variance ratio is fitted once under the null by
REML and each gene is then tested by generalised least squares — the FaST-LMM
approximation pyseer uses. Covariates: genome size and CDS count, because
susceptible strains carry measurably smaller genomes.

Significance is judged against a family-wise threshold from permuting the
phenotype, not Bonferroni: the features are heavily correlated, so the nominal
count of tests overstates the true multiplicity.

USAGE
-----
    python3 association.py --host-range /path/to/Host-Range --work work [--permutations 1000]
"""
import argparse, csv, os, sys, time, json
import numpy as np
from scipy import stats
from scipy.optimize import minimize_scalar


# ---------------------------------------------------------------- inputs ----

def assembly_metrics(host_range, strains, cache):
    if os.path.exists(cache):
        rows = list(csv.DictReader(open(cache), delimiter="\t"))
        by = {r["sample"]: r for r in rows}
        return np.array([[float(by[s]["contigs"]), float(by[s]["total_bp"]),
                          float(by[s]["n50"]), float(by[s]["cds"])] for s in strains])
    out = []
    for s in strains:
        d = os.path.join(host_range, "data/hosts/prokka_annotations", s)
        lens, cur = [], 0
        for line in open(os.path.join(d, f"{s}.fna")):
            if line.startswith(">"):
                if cur:
                    lens.append(cur)
                cur = 0
            else:
                cur += len(line.strip())
        if cur:
            lens.append(cur)
        lens.sort(reverse=True)
        total, acc, n50 = sum(lens), 0, 0
        for L in lens:
            acc += L
            if acc >= total / 2:
                n50 = L
                break
        cds = sum(1 for l in open(os.path.join(d, f"{s}.faa")) if l.startswith(">"))
        out.append(dict(sample=s, contigs=len(lens), total_bp=total, n50=n50, cds=cds))
    with open(cache, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", "contigs", "total_bp", "n50", "cds"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(out)
    return np.array([[r["contigs"], r["total_bp"], r["n50"], r["cds"]] for r in out], float)


def core_genome(host_range, strains, kin_cache, dist_cache):
    """Kinship and pairwise SNP distance from the core gene alignment."""
    if os.path.exists(kin_cache) and os.path.exists(dist_cache):
        return np.load(kin_cache), np.load(dist_cache)
    path = os.path.join(host_range,
                        "data/hosts/panaroo_pangenome/pangenome/core_gene_alignment.aln")
    names, seqs, cur = [], [], []
    for line in open(path):
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
                cur = []
            names.append(line[1:].strip())
        else:
            cur.append(line.strip())
    seqs.append("".join(cur))
    A = np.frombuffer("".join(seqs).encode(), dtype=np.uint8).reshape(len(names), -1)
    keep = np.isin(A, [ord(c) for c in "acgtACGT"]).all(0)
    A = A[:, keep]
    A = A[:, (A != A[0]).any(0)]
    order = [names.index(s) for s in strains]
    A = A[order]
    D = np.zeros((len(strains),) * 2)
    for i in range(len(strains)):
        D[i] = (A != A[i]).sum(1)
    major = np.array([np.bincount(A[:, j]).argmax() for j in range(A.shape[1])], dtype=A.dtype)
    Z = (A != major).astype(float)
    Z -= Z.mean(0)
    sd = Z.std(0)
    Z = Z[:, sd > 0] / sd[sd > 0]
    K = Z @ Z.T / Z.shape[1]
    np.save(kin_cache, K)
    np.save(dist_cache, D)
    return K, D


def load_recalled(path, order_path, strains):
    """Read either caller's output.

    recall_presence_absence.py:  cluster, query_len, bits      (no header)
    recall_genes.py:             family, n_clusters, query_len, bits  (header)
    """
    order = open(order_path).read().split("\n")
    idx = [order.index(s) for s in strains]
    names, mat, qlen = [], [], []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if f[0] == "family":
            continue
        if len(f) == 4:
            name, L, bits = f[0], int(f[2]), f[3]
        else:
            name, L, bits = f[0], int(f[1]), f[2]
        names.append(name)
        qlen.append(L)
        mat.append([int(x) for x in bits])
    return names, np.array(mat, dtype=np.int8)[:, idx], np.array(qlen)


# ----------------------------------------------------------------- model ----

def reml_delta(y, X, w, U):
    Uty, UtX = U.T @ y, U.T @ X
    def neg(logd):
        d = np.exp(logd)
        S = w + d
        A = UtX.T @ (UtX / S[:, None])
        try:
            b = np.linalg.solve(A, UtX.T @ (Uty / S))
        except np.linalg.LinAlgError:
            return 1e10
        r = Uty - UtX @ b
        n, p = len(y), X.shape[1]
        s = float((r ** 2 / S).sum() / (n - p))
        _, logdet = np.linalg.slogdet(A)
        return 0.5 * ((n - p) * np.log(s) + np.log(S).sum() + logdet)
    return float(np.exp(minimize_scalar(neg, bounds=(-10, 10), method="bounded").x))


def gls(y, X, g, w, U, delta):
    Xf = np.column_stack([X, g])
    S = w + delta
    UtX, Uty = U.T @ Xf, U.T @ y
    A = UtX.T @ (UtX / S[:, None])
    try:
        beta = np.linalg.solve(A, UtX.T @ (Uty / S))
        var = float((((Uty - UtX @ beta) ** 2 / S).sum() / (len(y) - Xf.shape[1]))
                    * np.linalg.inv(A)[-1, -1])
    except np.linalg.LinAlgError:
        return np.nan, np.nan
    if var <= 0:
        return np.nan, np.nan
    t = beta[-1] / np.sqrt(var)
    return float(beta[-1]), float(2 * stats.t.sf(abs(t), len(y) - Xf.shape[1]))


def scan(y, K, cov, M, names, min_count=5):
    w, U = np.linalg.eigh(K)
    w = np.maximum(w, 1e-9)
    X = np.column_stack([np.ones(len(y))] + ([cov] if cov is not None else []))
    delta = reml_delta(y, X, w, U)
    res = []
    for i in range(M.shape[0]):
        g = M[i].astype(float)
        c = g.sum()
        if c < min_count or c > len(y) - min_count:
            continue
        b, p = gls(y, X, g, w, U, delta)
        if np.isfinite(p):
            res.append((names[i], int(c), b, p))
    return delta, res


def lambda_gc(ps):
    ps = np.asarray([p for p in ps if np.isfinite(p) and p > 0])
    return float(np.median(stats.chi2.isf(ps, 1)) / stats.chi2.ppf(0.5, 1))


def bh(ps):
    ps = np.asarray(ps)
    o = np.argsort(ps)
    q = np.minimum.accumulate((ps[o] * len(ps) / np.arange(1, len(ps) + 1))[::-1])[::-1]
    out = np.empty(len(ps))
    out[o] = np.minimum(q, 1)
    return out


def permutation_threshold(y, K, cov, M, names, n_perm, seed=0):
    """Family-wise threshold: distribution of the smallest p under permuted phenotypes."""
    rng = np.random.default_rng(seed)
    mins = []
    for _ in range(n_perm):
        _, r = scan(rng.permutation(y), K, cov, M, names)
        mins.append(min(x[3] for x in r))
    return np.sort(mins)


# ------------------------------------------------------------------ main ----

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--host-range", required=True)
    ap.add_argument("--work", default="work")
    ap.add_argument("--calls", default="genes.tsv",
                    help="presence/absence calls inside --work")
    ap.add_argument("--permutations", type=int, default=300)
    ap.add_argument("--snp-threshold", type=int, default=10,
                    help="core SNPs within which strains count as one clone")
    args = ap.parse_args()
    W = args.work
    os.makedirs(W, exist_ok=True)

    pheno = list(csv.DictReader(
        open(os.path.join(args.host_range, "config/meta/host-range-results.tsv")),
        delimiter="\t"))
    strains = [r["sample"] for r in pheno]
    y = np.array([int(r["infection"]) for r in pheno], float)
    panel = np.array([s.split("-")[0] if s.split("-")[0] in ("BOR", "DTD", "VISA")
                      else "IN-HOUSE" for s in strains])

    asm = assembly_metrics(args.host_range, strains, os.path.join(W, "assembly_metrics.tsv"))
    K, D = core_genome(args.host_range, strains,
                       os.path.join(W, "kinship.npy"), os.path.join(W, "snp_distance.npy"))
    names, M, qlen = load_recalled(os.path.join(W, args.calls),
                                   os.path.join(W, "recall_strain_order.txt"), strains)
    z = lambda v: (v - v.mean()) / v.std()
    cov = np.column_stack([z(asm[:, 1]), z(asm[:, 3])])

    print(f"{len(strains)} strains, {int(y.sum())} susceptible | {len(names)} re-called clusters")
    n = len(strains)
    freq = M.sum(1)
    print(f"  now core (>=95%): {(freq >= 0.95 * n).sum()} | still testable: "
          f"{((freq >= 5) & (freq <= n - 5)).sum()}")

    delta, res = scan(y, K, cov, M, names)
    ps = np.array([r[3] for r in res])
    q = bh(ps)
    print(f"\nGWAS on re-called calls: tested={len(res)} lambda={lambda_gc(ps):.2f} delta={delta:.2f}")
    t0 = time.time()
    mins = permutation_threshold(y, K, cov, M, names, args.permutations)
    thr = mins[int(0.05 * len(mins))]
    print(f"{args.permutations} permutations in {time.time()-t0:.0f}s -> "
          f"family-wise 5% threshold p = {thr:.2e}")

    with open(os.path.join(W, "gwas.tsv"), "w") as fh:
        fh.write("cluster\tn_present\tbeta\tp\tq_bh\tfwer_p\n")
        for (g, c, b, p), qq in zip(res, q):
            fh.write(f"{g}\t{c}\t{b:.4f}\t{p:.3e}\t{qq:.4f}\t{(mins <= p).mean():.4f}\n")
    o = np.argsort(ps)
    print(f"\n{'cluster':40s}{'n+':>5s}{'beta':>8s}{'p':>10s}{'FWER':>8s}")
    for i in o[:10]:
        g, c, b, p = res[i]
        print(f"{g[:40]:40s}{c:5d}{b:8.3f}{p:10.2e}{(mins <= p).mean():8.3f}")
    print(f"significant at family-wise 5%: {(ps <= thr).sum()}")

    # ---- within-clone contrasts -------------------------------------------
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i in range(n):
        for j in range(i + 1, n):
            if D[i, j] <= args.snp_threshold:
                a, b_ = find(i), find(j)
                if a != b_:
                    parent[b_] = a
    clones = {}
    for i in range(n):
        clones.setdefault(find(i), []).append(i)
    mixed = {r: idx for r, idx in clones.items()
             if len(idx) > 1 and 0 < y[idx].sum() < len(idx)}
    print(f"\nclones at <= {args.snp_threshold} core SNPs: {len(clones)} groups; "
          f"{len(mixed)} contain both susceptible and resistant strains")
    with open(os.path.join(W, "within_clone.tsv"), "w") as fh:
        fh.write("clone\tn\tn_susceptible\tcluster\tpresent_in_susceptible\t"
                 "present_in_resistant\tp_fisher\n")
        rows = []
        for ci, (root, idx) in enumerate(sorted(mixed.items(), key=lambda kv: -len(kv[1])), 1):
            idx = np.array(idx)
            ys = y[idx]
            for gi_, gname in enumerate(names):
                g = M[gi_, idx]
                if g.sum() in (0, len(idx)):
                    continue
                a = int(((g == 1) & (ys == 1)).sum()); b_ = int(((g == 1) & (ys == 0)).sum())
                c_ = int(((g == 0) & (ys == 1)).sum()); d_ = int(((g == 0) & (ys == 0)).sum())
                p = stats.fisher_exact([[a, b_], [c_, d_]])[1]
                rows.append((ci, len(idx), int(ys.sum()), gname, f"{a}/{int(ys.sum())}",
                             f"{b_}/{int(len(idx)-ys.sum())}", p))
        rows.sort(key=lambda r: r[6])
        for r in rows:
            fh.write("\t".join(str(x) for x in r[:6]) + f"\t{r[6]:.4g}\n")
    print(f"  within-clone discriminating features written: {len(rows)}")
    if rows:
        print(f"  strongest: {rows[0][3]} p={rows[0][6]:.3g} "
              f"(present in {rows[0][4]} susceptible vs {rows[0][5]} resistant)")


if __name__ == "__main__":
    main()
