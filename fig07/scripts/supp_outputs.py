#!/usr/bin/env python3
"""Supplementary tables and the cluster-splitting demonstration for Figure 7.

Produces three things, all written into --work:

  supp_table_association.tsv
      Every gene family tested, with counts, effect, uncorrected p, BH q, the
      permutation family-wise p, and the p when the VISA collection is excluded.
      The last column is what shows that the nominally significant hits are
      confined to that collection.

  supp_table_alleles.tsv
      The WTA/LTA allele analysis, summarised per gene rather than per position:
      how many variants were tested, the strongest one, and whether anything
      cleared the family-wise threshold.

  cluster_splitting.tsv / cluster_splitting_stats.tsv
      The evidence that Panaroo's presence/absence matrix cannot support gene
      association claims here. For genes that Panaroo split across clusters, this
      gives each arm's presence pattern and each arm's association, next to the
      merged gene's. rplI is the case to look at: two arms of 24 and 59 strains,
      significant in OPPOSITE directions, for a ribosomal protein that is in fact
      present in 82 of 83 strains.

USAGE
    python3 supp_outputs.py --host-range /path/to/Host-Range --work work
"""
import argparse, csv, os, sys, collections
import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from association import (core_genome, assembly_metrics, load_recalled, scan, bh)

# Genes Panaroo split, with the arms named as they appear in the raw matrix.
SPLIT_CASES = {
    "rplI (50S ribosomal protein L9)": ["rplI", "group_273"],
    "sraP (serine-rich adhesin)": ["sraP~~~sraP_1~~~sraP_4",
                                   "sraP_1~~~sraP_2~~~sraP_3~~~sraP"],
    "agrB (accessory gene regulator B)": ["agrB", "agrB~~~agrB_1~~~agrB_2"],
    "metXA (homoserine O-acetyltransferase)": ["metXA~~~metXA_1",
                                               "metXA_2~~~metXA_1~~~metXA"],
}


def raw_matrix(host_range, strains):
    pan = os.path.join(host_range, "data/hosts/panaroo_pangenome/pangenome")
    rows = [l.rstrip("\n").split("\t") for l in
            open(os.path.join(pan, "gene_presence_absence.Rtab"))]
    order = rows[0][1:]
    idx = [order.index(s) for s in strains]
    return {r[0]: np.array([int(x) for x in r[1:]], float)[idx] for r in rows[1:]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--host-range", required=True)
    ap.add_argument("--work", default="work")
    ap.add_argument("--calls", default="genes.tsv")
    args = ap.parse_args()
    W = args.work

    pheno = list(csv.DictReader(
        open(os.path.join(args.host_range, "config/meta/host-range-results.tsv")),
        delimiter="\t"))
    strains = [r["sample"] for r in pheno]
    y = np.array([int(r["infection"]) for r in pheno], float)
    panel = np.array([s.split("-")[0] if s.split("-")[0] in ("BOR", "DTD", "VISA")
                      else "IN-HOUSE" for s in strains])
    K, D = core_genome(args.host_range, strains,
                       os.path.join(W, "kinship.npy"), os.path.join(W, "snp_distance.npy"))
    asm = assembly_metrics(args.host_range, strains, os.path.join(W, "assembly_metrics.tsv"))
    z = lambda v: (v - v.mean()) / v.std()
    cov = np.column_stack([z(asm[:, 1]), z(asm[:, 3])])
    names, M, _ = load_recalled(os.path.join(W, args.calls),
                                os.path.join(W, "recall_strain_order.txt"), strains)

    # ---- supplementary table 1: association results ------------------------
    _, res = scan(y, K, cov, M, names)
    ps = np.array([r[3] for r in res])
    q = bh(ps)
    keep = panel != "VISA"
    sel = np.where(keep)[0]
    _, res_nv = scan(y[sel], K[np.ix_(sel, sel)],
                     np.column_stack([z(asm[sel, 1]), z(asm[sel, 3])]),
                     M[:, sel], names)
    p_nv = {r[0]: r[3] for r in res_nv}

    mins_path = os.path.join(W, "permutation_minima.npy")
    mins = np.load(mins_path) if os.path.exists(mins_path) else None

    gi = {g: i for i, g in enumerate(names)}
    rows = []
    for (g, c, b, p), qq in zip(res, q):
        v = M[gi[g]]
        a = int(((v == 1) & (y == 0)).sum()); bb = int(((v == 1) & (y == 1)).sum())
        cc = int(((v == 0) & (y == 0)).sum()); dd = int(((v == 0) & (y == 1)).sum())
        orr = stats.fisher_exact([[a, bb], [cc, dd]])[0]
        rows.append(dict(
            gene_family=g, strains_with_gene=c,
            susceptible_with_gene=bb, susceptible_without_gene=dd,
            odds_ratio=round(orr, 3) if np.isfinite(orr) else "NA",
            beta=round(b, 4), p=f"{p:.3e}", q_bh=round(qq, 4),
            fwer_p=(f"{(mins <= p).mean():.4f}" if mins is not None else "NA"),
            p_excluding_VISA=(f"{p_nv[g]:.3e}" if g in p_nv else "NA")))
    rows.sort(key=lambda r: float(r["p"]))
    out1 = os.path.join(W, "supp_table_association.tsv")
    with open(out1, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print(f"wrote {out1}  ({len(rows)} gene families)")

    # ---- supplementary table 2: alleles, summarised per gene ---------------
    allele_path = os.path.join(W, "allele_association.tsv")
    if os.path.exists(allele_path):
        by_gene = collections.defaultdict(list)
        for r in csv.DictReader(open(allele_path), delimiter="\t"):
            by_gene[r["feature"].split(":")[0]].append(r)
        out2 = os.path.join(W, "supp_table_alleles.tsv")
        with open(out2, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["gene", "features_tested", "strongest_feature",
                        "beta", "p", "fwer_p", "significant_family_wise"])
            for g, rs in sorted(by_gene.items()):
                best = min(rs, key=lambda r: float(r["p"]))
                w.writerow([g, len(rs), best["feature"], best["beta"], best["p"],
                            best["fwer_p"],
                            "yes" if float(best["fwer_p"]) < 0.05 else "no"])
        print(f"wrote {out2}  ({len(by_gene)} genes)")
    else:
        print(f"no {allele_path}; skipping the allele table")

    # ---- cluster splitting: the evidence the raw matrix is unusable --------
    raw = raw_matrix(args.host_range, strains)
    per_strain, stat_rows = [], []
    for label, arms in SPLIT_CASES.items():
        arms_present = [a for a in arms if a in raw]
        if len(arms_present) < 2:
            print(f"  {label}: arms not found in the raw matrix, skipped")
            continue
        merged = np.clip(sum(raw[a] for a in arms_present), 0, 1)
        feats = [raw[a] for a in arms_present] + [merged]
        labels = [f"{label}|arm {i+1}" for i in range(len(arms_present))] + [f"{label}|merged"]
        _, rs = scan(y, K, cov, np.array(feats), labels, min_count=0)
        for (lab, c, b, p), vec in zip(rs, feats):
            stat_rows.append(dict(
                gene=label, feature=lab.split("|")[1],
                cluster=(arms_present + ["merged"])[labels.index(lab)],
                strains_with=int(vec.sum()),
                susceptible_with=int(((vec == 1) & (y == 1)).sum()),
                beta=round(b, 4), p=f"{p:.3e}",
                direction=("present in susceptible" if b > 0 else "present in resistant")))
        for si, s in enumerate(strains):
            row = dict(strain=s, gene=label, susceptible=int(y[si]), panel=panel[si],
                       merged=int(merged[si]))
            for ai, a in enumerate(arms_present):
                row[f"arm{ai+1}"] = int(raw[a][si])
            per_strain.append(row)

    if stat_rows:
        out3 = os.path.join(W, "cluster_splitting_stats.tsv")
        with open(out3, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(stat_rows[0]), delimiter="\t")
            w.writeheader(); w.writerows(stat_rows)
        out4 = os.path.join(W, "cluster_splitting.tsv")
        fields = ["strain", "gene", "susceptible", "panel", "arm1", "arm2", "merged"]
        with open(out4, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            w.writeheader(); w.writerows(per_strain)
        print(f"wrote {out3} and {out4}")
        print(f"\n{'gene':40s}{'feature':10s}{'n':>4s}{'beta':>9s}{'p':>11s}  direction")
        for r in stat_rows:
            print(f"{r['gene'][:40]:40s}{r['feature']:10s}{r['strains_with']:4d}"
                  f"{r['beta']:9.3f}{r['p']:>11s}  {r['direction']}")


if __name__ == "__main__":
    main()
