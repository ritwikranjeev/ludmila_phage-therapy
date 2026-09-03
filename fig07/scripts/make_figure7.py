#!/usr/bin/env python3
"""Figure 7 (host range) and its supplementary panel (detection limit).

MAIN FIGURE — descriptive, no association claim.
  7a  Neighbour-joining tree of the 83 hosts from core-genome SNP distances,
      tips marked with susceptibility group and source collection. Supports the
      statement that susceptibility is not lineage-associated: relatedness
      explains 7% of the variance and susceptible and resistant strains are
      interleaved across the tree.
  7b  Susceptibility by source collection, which is where the structure actually
      sits (VISA 21% susceptible against ~70% elsewhere).

SUPPLEMENTARY — what the panel could have detected.
  Power to detect a gene-content association at the permutation-derived
  threshold, by gene frequency and odds ratio. This is what makes the negative
  result quantitative rather than an absence of evidence.

USAGE
    python3 make_figure7.py --host-range /path/to/Host-Range --work work --out output
"""
import argparse, csv, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from association import core_genome, assembly_metrics, reml_delta, gls

GROUP_COLOUR = {"A": "#2c7fb8", "B": "#f0a63a", "C": "#d1495b"}
GROUP_LABEL = {"A": "A  fully susceptible", "B": "B  partially resistant",
               "C": "C  fully resistant"}
PANEL_MARK = {"BOR": "o", "DTD": "s", "VISA": "^", "IN-HOUSE": "D"}


# ------------------------------------------------------------------ tree ----

def neighbour_joining(D, labels):
    """Saitou-Nei neighbour joining. Returns (children, heights) for plotting."""
    D = D.astype(float).copy()
    n = len(labels)
    nodes = list(range(n))
    children = {}
    lengths = {i: 0.0 for i in range(n)}
    next_id = n
    active = list(range(n))
    dist = {(i, j): D[i, j] for i in range(n) for j in range(n)}

    def d(a, b):
        return dist[(a, b)] if (a, b) in dist else dist[(b, a)]

    while len(active) > 2:
        m = len(active)
        r = {i: sum(d(i, k) for k in active if k != i) / (m - 2) for i in active}
        best, bi, bj = None, None, None
        for ii in range(m):
            for jj in range(ii + 1, m):
                i, j = active[ii], active[jj]
                q = d(i, j) - r[i] - r[j]
                if best is None or q < best:
                    best, bi, bj = q, i, j
        li = 0.5 * d(bi, bj) + 0.5 * (r[bi] - r[bj])
        lj = d(bi, bj) - li
        children[next_id] = [(bi, max(li, 0.0)), (bj, max(lj, 0.0))]
        for k in active:
            if k in (bi, bj):
                continue
            dist[(next_id, k)] = 0.5 * (d(bi, k) + d(bj, k) - d(bi, bj))
            dist[(k, next_id)] = dist[(next_id, k)]
        active = [k for k in active if k not in (bi, bj)] + [next_id]
        next_id += 1
    i, j = active
    children[next_id] = [(i, max(d(i, j), 0.0)), (j, 0.0)]
    return children, next_id


def layout(children, root, n_tips):
    """Assign x (distance from root) and y (tip order) to every node."""
    x, y = {}, {}
    order = []

    def walk(node, depth):
        x[node] = depth
        if node not in children:
            order.append(node)
            y[node] = len(order) - 1
            return
        for child, length in children[node]:
            walk(child, depth + length)
        y[node] = np.mean([y[c] for c, _ in children[node]])

    walk(root, 0.0)
    return x, y, order


def draw_tree(ax, children, root, x, y, order, labels, groups, panels):
    for node, kids in children.items():
        ys = [y[c] for c, _ in kids]
        ax.plot([x[node]] * 2, [min(ys), max(ys)], color="0.35", lw=0.7, zorder=1)
        for c, _ in kids:
            ax.plot([x[node], x[c]], [y[c], y[c]], color="0.35", lw=0.7, zorder=1)
    xmax = max(x[t] for t in order)
    for t in order:
        g = groups[t]
        ax.scatter(x[t] + 0.012 * xmax, y[t], s=26, color=GROUP_COLOUR[g],
                   marker=PANEL_MARK[panels[t]], edgecolor="0.2", linewidth=0.4, zorder=3)
        ax.text(x[t] + 0.03 * xmax, y[t], labels[t], va="center", fontsize=5.2,
                color="0.15")
    ax.set_xlim(-0.01 * xmax, xmax * 1.22)
    ax.set_ylim(-1, len(order))
    ax.axis("off")


# ------------------------------------------------------------------ power ---

def power_curve(host_range, strains, K, asm, cache, threshold, reps, seed=0):
    if os.path.exists(cache):
        rows = list(csv.DictReader(open(cache), delimiter="\t"))
        return rows
    z = lambda v: (v - v.mean()) / v.std()
    X = np.column_stack([np.ones(len(strains)), z(asm[:, 1]), z(asm[:, 3])])
    w, U = np.linalg.eigh(K)
    w = np.maximum(w, 1e-9)
    rng = np.random.default_rng(seed)
    base = 53 / 83
    out = []
    for freq in (0.15, 0.30, 0.50):
        for OR in (2, 3, 5, 8, 12, 20, 40):
            hits = 0
            for _ in range(reps):
                g = (rng.random(len(strains)) < freq).astype(float)
                if not (5 <= g.sum() <= len(strains) - 5):
                    continue
                odds = base / (1 - base) * np.where(g == 1, 1 / OR, 1.0)
                y = (rng.random(len(strains)) < odds / (1 + odds)).astype(float)
                if not (5 <= y.sum() <= len(strains) - 5):
                    continue
                delta = reml_delta(y, X, w, U)
                _, p = gls(y, X, g, w, U, delta)
                hits += int(np.isfinite(p) and p <= threshold)
            out.append(dict(frequency=freq, odds_ratio=OR, power=hits / reps))
            print(f"  freq {freq:.2f}  OR {OR:2d}  power {hits/reps:.2f}", flush=True)
    with open(cache, "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=["frequency", "odds_ratio", "power"],
                             delimiter="\t")
        wtr.writeheader()
        wtr.writerows(out)
    return out


# ------------------------------------------------------------------- main ---

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--host-range", required=True)
    ap.add_argument("--work", default="work")
    ap.add_argument("--out", default="output")
    ap.add_argument("--threshold", type=float, default=1.84e-4,
                    help="family-wise 5%% threshold from association.py")
    ap.add_argument("--power-reps", type=int, default=150)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    pheno = list(csv.DictReader(
        open(os.path.join(args.host_range, "config/meta/host-range-results.tsv")),
        delimiter="\t"))
    strains = [r["sample"] for r in pheno]
    groups = {i: pheno[i]["group"] for i in range(len(strains))}
    panels = {i: (strains[i].split("-")[0] if strains[i].split("-")[0] in PANEL_MARK
                  else "IN-HOUSE") for i in range(len(strains))}
    K, D = core_genome(args.host_range, strains,
                       os.path.join(args.work, "kinship.npy"),
                       os.path.join(args.work, "snp_distance.npy"))
    asm = assembly_metrics(args.host_range, strains,
                           os.path.join(args.work, "assembly_metrics.tsv"))

    print("building neighbour-joining tree", flush=True)
    children, root = neighbour_joining(D, strains)
    x, y, order = layout(children, root, len(strains))

    # ---- Figure 7 ----------------------------------------------------------
    fig = plt.figure(figsize=(9.5, 10))
    gs = fig.add_gridspec(2, 1, height_ratios=[4.4, 1.0], hspace=0.16)
    ax = fig.add_subplot(gs[0])
    draw_tree(ax, children, root, x, y, order, strains, groups, panels)
    ax.set_title("Core-genome relatedness and OCPT-Ludmila susceptibility",
                 fontsize=11, weight="bold", loc="left")
    handles = [Line2D([], [], marker="o", ls="", color=GROUP_COLOUR[g],
                      markeredgecolor="0.2", markersize=6, label=GROUP_LABEL[g])
               for g in "ABC"]
    handles += [Line2D([], [], marker=PANEL_MARK[p], ls="", color="0.75",
                       markeredgecolor="0.2", markersize=6, label=p)
                for p in ["BOR", "DTD", "VISA", "IN-HOUSE"]]
    ax.legend(handles=handles, loc="lower left", frameon=False, fontsize=7,
              ncol=2, handletextpad=0.4)

    ax2 = fig.add_subplot(gs[1])
    panel_names = ["BOR", "DTD", "VISA", "IN-HOUSE"]
    counts = {p: [sum(1 for i in range(len(strains))
                      if panels[i] == p and groups[i] == g) for g in "ABC"]
              for p in panel_names}
    bottom = np.zeros(len(panel_names))
    for gi, g in enumerate("ABC"):
        vals = np.array([counts[p][gi] for p in panel_names], float)
        ax2.bar(panel_names, vals, bottom=bottom, color=GROUP_COLOUR[g],
                edgecolor="white", linewidth=0.6, label=GROUP_LABEL[g])
        bottom += vals
    for i, p in enumerate(panel_names):
        tot = sum(counts[p])
        ax2.text(i, tot + 0.8, f"{100*counts[p][0]/tot:.0f}% group A\n(n={tot})",
                 ha="center", fontsize=7, color="0.25")
    ax2.set_ylabel("strains", fontsize=8)
    ax2.set_ylim(0, max(sum(counts[p]) for p in panel_names) * 1.35)
    ax2.tick_params(labelsize=8)
    for side in ("top", "right"):
        ax2.spines[side].set_visible(False)
    ax2.set_title("Susceptibility by source collection", fontsize=10,
                  weight="bold", loc="left")

    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(args.out, f"figure7_host_range.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote figure7_host_range.png/.pdf", flush=True)

    # ---- Supplementary: detection limit ------------------------------------
    print("simulating power", flush=True)
    rows = power_curve(args.host_range, strains, K, asm,
                       os.path.join(args.work, "power_curve.tsv"),
                       args.threshold, args.power_reps)
    fig2, axp = plt.subplots(figsize=(6.2, 4.2))
    for freq, marker in zip(("0.15", "0.3", "0.5"), ("o", "s", "^")):
        sel = [r for r in rows if f"{float(r['frequency']):.2f}" == f"{float(freq):.2f}"]
        sel.sort(key=lambda r: float(r["odds_ratio"]))
        axp.plot([float(r["odds_ratio"]) for r in sel],
                 [float(r["power"]) for r in sel], marker=marker, lw=1.4,
                 label=f"gene in {int(round(float(freq)*len(strains)))} of {len(strains)} strains")
    axp.axhline(0.8, ls="--", lw=0.8, color="0.5")
    axp.text(2.1, 0.82, "80% power", fontsize=7, color="0.4")
    axp.set_xscale("log")
    axp.set_xticks([2, 3, 5, 8, 12, 20, 40])
    axp.set_xticklabels(["2", "3", "5", "8", "12", "20", "40"])
    axp.set_xlabel("odds ratio", fontsize=9)
    axp.set_ylabel("power at the family-wise 5% threshold", fontsize=9)
    axp.set_ylim(0, 1.02)
    axp.legend(frameon=False, fontsize=8, loc="upper left")
    for side in ("top", "right"):
        axp.spines[side].set_visible(False)
    axp.set_title("Effect sizes this panel could detect", fontsize=10,
                  weight="bold", loc="left")
    for ext in ("png", "pdf"):
        fig2.savefig(os.path.join(args.out, f"supp_figure7_detection_limit.{ext}"),
                     dpi=300, bbox_inches="tight")
    plt.close(fig2)
    print("wrote supp_figure7_detection_limit.png/.pdf", flush=True)


if __name__ == "__main__":
    main()
