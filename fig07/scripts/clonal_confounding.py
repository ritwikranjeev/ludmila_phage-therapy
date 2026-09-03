#!/usr/bin/env python3
"""Why the three genome-wide-significant families are not a finding.

WHAT THE SCAN REPORTED
----------------------
Three gene families clear both Benjamini-Hochberg and the permutation
family-wise threshold in the 83-strain panel:

    group_2749   OR 13.7   p 2.4e-05   q 0.031   fwer 0.006
    group_2693   OR  6.1   p 9.2e-05   q 0.046   fwer 0.020
    group_2736   OR  8.5   p 1.1e-04   q 0.046   fwer 0.028

Taken at face value that is a positive result. It is not, and this script is
the demonstration.

THE THREE STEPS
---------------
1. THEY ARE ONE LOCUS, NOT THREE FINDINGS. Pairwise Jaccard between the four
   top families is 0.77-0.94; they are co-inherited. Treating them as three
   independent hits triple-counts one observation.

2. THE LOCUS MARKS THE VISA SUBPANEL. It is present in 13 of 14 VISA strains.
   The VISA strains are themselves markedly phage-resistant (21% susceptible
   vs 72% for everything else, Fisher p = 5e-04), so any marker of that
   subpanel inherits the association.

3. THE VISA SUBPANEL IS TWO CLONES. At a 50-SNP threshold the 14 VISA isolates
   collapse to 2 lineages. The effective sample size behind the association is
   therefore closer to 2 than to 14, and the p-value computed over 83 nominally
   independent strains is correspondingly inflated.

Collapsing the WHOLE panel to independent lineages and testing once per lineage
removes the association at every threshold tried (p = 1.00, 0.62, 0.58).

WHY THE MIXED MODEL DID NOT CATCH IT
------------------------------------
Two reasons, and both are worth stating in the Methods. Core-genome kinship
downweights relatedness but does not eliminate it when the signal is carried by
a small number of large clones. And the permutation threshold shuffles the
phenotype freely across all 83 strains, which destroys the clonal correlation
the real data has -- so the null it builds is easier to beat than the truth,
making the family-wise threshold anti-conservative here.

PHENOTYPE ENCODING
------------------
Also checked, because the pipeline defines several: susceptibility as
discrete_scale == 3 (53/83, what the paper uses) versus discrete_scale > 0
(72/83). The locus association is p = 0.0033 and 0.0028 respectively. The
encoding does not matter.

USAGE
    python3 clonal_confounding.py --work work --host-range /path/to/Host-Range
"""
import argparse, csv, os
import numpy as np
from scipy import stats

TOP = ["group_2749", "group_2693", "group_2736", "group_2730"]


def connected_components(D, thr):
    """Single-linkage clustering of strains at a SNP threshold."""
    n = D.shape[0]
    lab = np.full(n, -1)
    g = 0
    for i in range(n):
        if lab[i] != -1:
            continue
        stack, seen = [i], set()
        while stack:
            k = stack.pop()
            if k in seen:
                continue
            seen.add(k)
            lab[k] = g
            for j in range(n):
                if j not in seen and D[k, j] <= thr:
                    stack.append(j)
        g += 1
    return lab, g


def fisher(a, b, c, d):
    return stats.fisher_exact([[a, b], [c, d]])[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default="work")
    ap.add_argument("--host-range", required=True)
    ap.add_argument("--thresholds", default="50,100,500")
    args = ap.parse_args()
    W = args.work

    strains = [l.strip() for l in open(os.path.join(W, "recall_strain_order.txt")) if l.strip()]
    D = np.load(os.path.join(W, "snp_distance.npy"))
    ph = {r["sample"]: r for r in csv.DictReader(
        open(os.path.join(args.host_range, "config/meta/host-range-results.tsv")), delimiter="\t")}
    y = np.array([int(ph[s]["infection"]) for s in strains])
    disc = np.array([int(ph[s]["discrete_scale"]) for s in strains])
    visa = np.array([s.startswith("VISA") for s in strains])

    calls = {r["family"]: np.array([int(c) for c in r["calls"]])
             for r in csv.DictReader(open(os.path.join(W, "genes.tsv")), delimiter="\t")}
    M = np.vstack([calls[t] for t in TOP])
    locus = (M.sum(0) == len(TOP)).astype(int)

    rows = []

    # 1. are the top families one locus?
    print("1. co-inheritance of the top families (Jaccard)")
    for i in range(len(TOP)):
        for j in range(i + 1, len(TOP)):
            a, b = calls[TOP[i]], calls[TOP[j]]
            jac = (a & b).sum() / (a | b).sum()
            print(f"     {TOP[i]:12s} vs {TOP[j]:12s}  {jac:.3f}")
            rows.append(dict(test="jaccard", detail=f"{TOP[i]}|{TOP[j]}", value=f"{jac:.4f}", p=""))

    # 2. does the locus mark the VISA subpanel, and is that subpanel resistant?
    print("\n2. the VISA subpanel")
    p_visa = fisher(int(y[visa].sum()), int((~y.astype(bool))[visa].sum()),
                    int(y[~visa].sum()), int((~y.astype(bool))[~visa].sum()))
    print(f"     locus present in {int(locus[visa].sum())}/{int(visa.sum())} VISA "
          f"and {int(locus[~visa].sum())}/{int((~visa).sum())} non-VISA")
    print(f"     VISA susceptible {int(y[visa].sum())}/{int(visa.sum())} vs "
          f"non-VISA {int(y[~visa].sum())}/{int((~visa).sum())};  Fisher p = {p_visa:.4g}")
    rows.append(dict(test="visa_resistance", detail="VISA vs non-VISA susceptibility",
                     value=f"{int(y[visa].sum())}/{int(visa.sum())} vs {int(y[~visa].sum())}/{int((~visa).sum())}",
                     p=f"{p_visa:.4g}"))

    # 3. how many independent lineages, and does the association survive?
    print("\n3. collapsing to independent lineages")
    for thr in [int(t) for t in args.thresholds.split(",")]:
        lab, g = connected_components(D, thr)
        ry, rl = [], []
        for c in range(g):
            m = lab == c
            ry.append(int(y[m].mean() >= 0.5))
            rl.append(int(locus[m].mean() >= 0.5))
        ry, rl = np.array(ry), np.array(rl)
        p = fisher(int(((rl == 1) & (ry == 1)).sum()), int(((rl == 1) & (ry == 0)).sum()),
                   int(((rl == 0) & (ry == 1)).sum()), int(((rl == 0) & (ry == 0)).sum()))
        nvisa = connected_components(D[np.ix_(visa, visa)], thr)[1]
        print(f"     <= {thr:4d} SNPs: {g:2d} lineages overall, {nvisa} among VISA"
              f"   locus association p = {p:.3f}")
        rows.append(dict(test="lineage_collapsed", detail=f"<= {thr} SNPs, {g} lineages",
                         value=f"{nvisa} VISA lineages", p=f"{p:.4g}"))

    # 4. phenotype encoding
    print("\n4. phenotype encoding")
    for name, pheno in (("discrete_scale == 3", y), ("discrete_scale > 0", (disc > 0).astype(int))):
        p = fisher(int((locus & pheno).sum()), int((locus & (1 - pheno)).sum()),
                   int(((1 - locus) & pheno).sum()), int(((1 - locus) & (1 - pheno)).sum()))
        print(f"     {name:22s} {int(pheno.sum())}/83 susceptible   locus p = {p:.4g}")
        rows.append(dict(test="phenotype_encoding", detail=name,
                         value=f"{int(pheno.sum())}/83 susceptible", p=f"{p:.4g}"))

    out = os.path.join(W, "clonal_confounding.tsv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["test", "detail", "value", "p"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
