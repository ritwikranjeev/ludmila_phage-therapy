#!/usr/bin/env python3
"""Build the Figure 7 gene-presence grid and its per-gene statistics.

WHAT THIS REPLACES
------------------
The original figure drew its gene grid from
`Host-Range/results/01_feature_tables/merged_features.tsv`, which holds raw
Panaroo cluster calls. Those calls are wrong for several genes in a way that is
visible to anyone who knows the pathway:

    gene    merged_features.tsv    re-called
    tagO             3/83             83/83
    tarF             3/83             83/83
    tarK             3/83             78/83
    tarM             3/83             33/83
    clfB            51/83             75/83

tagO is essential for wall teichoic acid synthesis; 3/83 is a calling failure,
not a finding. The same cluster-splitting that produced those numbers is what
manufactured the sraP association (see supp_figure_cluster_splitting.R). So the
grid here is drawn from the corrected gene families in work/genes.tsv instead,
mapped back onto the original candidate gene names.

WHAT IS SHOWN
-------------
The 61-gene candidate panel from the original figure covers wall teichoic acid
synthesis, LTA, the dlt operon, surface adhesins and capsule. Of those:

  * 33 map to a re-called family. The 28 that do not are almost entirely the
    capsule locus (capA-capO), which reads 0/83 in the original table as well,
    so it is either genuinely absent from this panel or was never called.
  * 30 of the 33 are variable. The 3 invariant ones are dropped: an all-present
    or all-absent column carries no information and reads as a solid bar.

Two caveats travel with the output rather than being silently resolved:

  * tarK and tarB merge into a single re-called family, so they would render as
    identical rows. Only the merged family is emitted, labelled as such.
  * ltaS reads 2/83 in BOTH call sets. LtaS is essential in S. aureus, so this
    is implausible either way and is not explained by cluster splitting. It is
    flagged in the output rather than dropped.

STATISTICS
----------
Each gene is tested with the same mixed model as the genome-wide scan --
core-genome kinship, genome size and CDS count as covariates -- so the numbers
in the figure are on the same footing as those in supp_table_association.tsv.
Benjamini-Hochberg is applied across the genes shown. An uncorrected Fisher
exact p is reported alongside for readers who want the raw contingency test.

USAGE
    python3 gene_grid.py --host-range /path/to/Host-Range --work work
"""
import argparse, csv, os, re, sys
import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from association import core_genome, assembly_metrics, scan, bh

PANEL = {
    "WTA synthesis": ["tarO", "tagO", "tarA", "tagA", "tarB", "tagB", "tarD", "tagD",
                      "tarF", "tagF", "tarL", "tarK", "tagH", "tagG", "tagX"],
    "WTA glycosylation": ["tarM", "tarS", "tarP", "tarN", "tarQ", "tarV", "gttA", "yfhO", "gtlB"],
    "LTA": ["ltaS", "ltaA", "pgcA", "gtaB", "ypfP", "ugtP"],
    "dlt operon": ["mprF", "dltA", "dltB", "dltC", "dltD"],
    "Surface adhesins": ["sdrC", "sdrD", "sdrE", "clfA", "clfB", "fnbA", "fnbB",
                         "spa", "sasA", "sasG", "sraP"],
    "Capsule": ["capA", "capB", "capC", "capD", "capE", "capF", "capG", "capH",
                "capI", "capJ", "capK", "capL", "capM", "capN", "capO"],
}
IMPLAUSIBLE = {"ltaS": "essential in S. aureus; 2/83 in both call sets, cause unresolved"}


def base_name(n):
    return re.split(r"[_~]", n)[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--host-range", required=True)
    ap.add_argument("--work", default="work")
    args = ap.parse_args()
    W = args.work

    strains = [l.strip() for l in open(os.path.join(W, "recall_strain_order.txt")) if l.strip()]
    ph = {r["sample"]: r for r in csv.DictReader(
        open(os.path.join(args.host_range, "config/meta/host-range-results.tsv")), delimiter="\t")}
    y = np.array([int(ph[s]["infection"]) for s in strains], float)

    calls = {r["family"]: np.array([int(c) for c in r["calls"]])
             for r in csv.DictReader(open(os.path.join(W, "genes.tsv")), delimiter="\t")}
    fams = {r["family"]: r["clusters"]
            for r in csv.DictReader(open(os.path.join(W, "gene_families.tsv")), delimiter="\t")}

    index = {}
    for f in calls:
        index.setdefault(base_name(f), []).append(f)
    for f, cl in fams.items():
        for c in cl.split(","):
            index.setdefault(base_name(c.strip()), []).append(f)

    # map candidate names onto re-called families, keeping the panel grouping
    chosen, seen_family = [], {}
    unmapped, invariant, merged = [], [], []
    for group, genes in PANEL.items():
        for g in genes:
            hits = sorted(set(index.get(g, [])))
            hits = [h for h in hits if h in calls]
            if not hits:
                unmapped.append(g)
                continue
            fam = hits[0]
            v = calls[fam]
            if v.sum() in (0, len(strains)):
                invariant.append(g)
                continue
            if fam in seen_family:
                merged.append((g, seen_family[fam]))
                continue
            seen_family[fam] = g
            chosen.append(dict(gene=g, group=group, family=fam, vec=v))

    print(f"candidate panel: {sum(len(v) for v in PANEL.values())} genes")
    print(f"  mapped and variable : {len(chosen)}")
    print(f"  invariant (dropped) : {len(invariant)}  {invariant}")
    print(f"  unmapped (dropped)  : {len(unmapped)}")
    if merged:
        print(f"  merged into another family: " +
              ", ".join(f"{a} -> shown as {b}" for a, b in merged))

    # ---- statistics, same model as the genome-wide scan ----------------------
    K, D = core_genome(args.host_range, strains,
                       os.path.join(W, "kinship.npy"), os.path.join(W, "snp_distance.npy"))
    asm = assembly_metrics(args.host_range, strains, os.path.join(W, "assembly_metrics.tsv"))
    z = lambda v: (v - v.mean()) / v.std()
    cov = np.column_stack([z(asm[:, 1]), z(asm[:, 3])])

    M = np.vstack([c["vec"] for c in chosen]).astype(float)
    names = [c["gene"] for c in chosen]
    delta, res = scan(y, K, cov, M, names, min_count=1)
    lmm = {name: (beta, p) for name, count, beta, p in res}
    print(f"  variance ratio delta = {delta:.4g}; {len(res)}/{len(names)} genes testable")

    rows = []
    for c in chosen:
        v = c["vec"]
        a = int((v & (y > 0)).sum()); b = int(v.sum()) - a
        cc = int(((1 - v) & (y > 0)).sum()); d = len(strains) - a - b - cc
        fisher = stats.fisher_exact([[a, b], [cc, d]])[1]
        beta, p = lmm.get(c["gene"], (np.nan, np.nan))
        rows.append(dict(gene=c["gene"], panel_group=c["group"], family=c["family"],
                         n_present=int(v.sum()), susceptible_with_gene=a,
                         susceptible_without_gene=cc,
                         beta=f"{beta:.4g}" if beta == beta else "NA",
                         p_lmm=f"{p:.4g}" if p == p else "NA",
                         p_fisher=f"{fisher:.4g}",
                         note=IMPLAUSIBLE.get(c["gene"], "")))

    ps = np.array([float(r["p_lmm"]) if r["p_lmm"] != "NA" else np.nan for r in rows])
    ok = ~np.isnan(ps)
    q = np.full(len(rows), np.nan)
    if ok.sum():
        q[ok] = bh(ps[ok])
    for r, qq in zip(rows, q):
        r["q_bh"] = f"{qq:.4g}" if qq == qq else "NA"

    out_stats = os.path.join(W, "gene_grid_stats.tsv")
    with open(out_stats, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    out_grid = os.path.join(W, "gene_grid.tsv")
    with open(out_grid, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["strain", "gene", "panel_group", "present"])
        for c in chosen:
            for s, v in zip(strains, c["vec"]):
                w.writerow([s, c["gene"], c["group"], int(v)])

    sig = [r for r in rows if r["q_bh"] != "NA" and float(r["q_bh"]) < 0.05]
    print(f"\nwrote {out_grid} and {out_stats}")
    print(f"genes with q < 0.05 after BH across the {len(rows)} shown: {len(sig)}")
    for r in sorted(rows, key=lambda r: float(r["p_lmm"]) if r["p_lmm"] != "NA" else 1)[:6]:
        print(f"  {r['gene']:7s} {r['n_present']:3d}/83  p={r['p_lmm']:>9s}  q={r['q_bh']:>8s}")


if __name__ == "__main__":
    main()
