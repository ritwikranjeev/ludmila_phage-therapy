#!/usr/bin/env python3
"""
Mauve-style pre vs post figure: synteny ribbons from nucmer, genome bars, and
labels for genes hit by variants (e.g. unique-to-post rows in mutations.tsv).
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def read_fasta_contig_length(fasta: Path, contig: str) -> int:
    name = None
    length = 0
    with fasta.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if name == contig:
                    return length
                name = line[1:].split()[0]
                length = 0
            elif name == contig:
                length += len(line.strip())
    if name == contig:
        return length
    raise SystemExit(f"Contig {contig!r} not found in {fasta}")


def parse_show_coords_rows(
    delta: Path, ref_contig: str, qry_contig: str
) -> list[tuple[int, int, int, int]]:
    sc = "show-coords"
    prefix = os.environ.get("CONDA_PREFIX", "")
    if prefix:
        cand = Path(prefix) / "bin" / "show-coords"
        if cand.is_file():
            sc = str(cand)
    proc = subprocess.run(
        [sc, "-H", "-c", "-l", str(delta)],
        check=True,
        capture_output=True,
        text=True,
    )
    rows: list[tuple[int, int, int, int]] = []
    for line in proc.stdout.splitlines():
        line = line.strip()
        if not line or line[0] not in "0123456789":
            continue
        parts = [x.strip() for x in line.split("|")]
        if len(parts) < 7:
            continue
        tail = parts[-1].split()
        if len(tail) < 2:
            continue
        ref_id, qry_id = tail[-2], tail[-1]
        if ref_id != ref_contig or qry_id != qry_contig:
            continue
        try:
            r1, r2 = (int(x) for x in parts[0].split())
            q1, q2 = (int(x) for x in parts[1].split())
        except ValueError:
            continue
        rows.append((r1, r2, q1, q2))
    rows.sort(key=lambda t: (min(t[0], t[1]), max(t[0], t[1])))
    return rows


def ref_to_qry(
    pos: int, segments: list[tuple[int, int, int, int]]
) -> int | None:
    for r1, r2, q1, q2 in segments:
        lo_r, hi_r = min(r1, r2), max(r1, r2)
        if not (lo_r <= pos <= hi_r):
            continue
        if hi_r == lo_r:
            return int((q1 + q2) / 2)
        t = (pos - r1) / (r2 - r1)
        return int(round(q1 + t * (q2 - q1)))
    return None


def parse_gff_genes(
    gff: Path, seq: str, types: frozenset[str]
) -> dict[str, tuple[int, int, str]]:
    """locus_tag -> (start, end, strand) on seq; 1-based inclusive."""
    out: dict[str, tuple[int, int, str]] = {}
    tag_re = re.compile(r"locus_tag=([^;\s]+)")
    with gff.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] != seq or f[2] not in types:
                continue
            start, end = int(f[3]), int(f[4])
            strand = f[6] if f[6] in "+-" else "+"
            m = tag_re.search(f[8])
            if not m:
                continue
            out[m.group(1)] = (start, end, strand)
    return out


def load_mutation_genes(
    path: Path,
    chrom: str,
    only_unique_post: bool,
) -> dict[str, list[tuple[int, str]]]:
    """
    gene -> [(pos, short_note), ...] from TSV.
    """
    by_gene: dict[str, list[tuple[int, str]]] = defaultdict(list)
    with path.open(newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            if row.get("chrom") != chrom:
                continue
            if only_unique_post:
                u = (row.get("unique_to_post") or "").strip().lower()
                if u not in ("true", "1", "yes"):
                    continue
            gene = (row.get("gene") or "").strip()
            if not gene or gene == ".":
                continue
            try:
                pos = int(row["pos"])
            except (KeyError, ValueError):
                continue
            eff = (row.get("effect") or "").strip()
            detail = (row.get("effect_detail") or "").strip()
            vtype = (row.get("variant_type") or "").strip()
            note = f"{vtype}" + (f", {detail}" if detail else "")
            if eff:
                note = f"{eff}: {note}"
            by_gene[gene].append((pos, note))
    return dict(by_gene)


def gene_arrow_patches(
    start: int,
    end: int,
    strand: str,
    y: float,
    height: float,
    x0: float,
    scale: float,
    facecolor: str,
    edgecolor: str,
    lw: float,
) -> list[mpatches.Patch]:
    left = x0 + (start - 1) * scale
    right = x0 + end * scale
    span = max(right - left, scale * 150)
    body_h = height * 0.55
    head_w = min(span * 0.14, scale * 9000)
    head_w = max(head_w, scale * 350)
    head_w = min(head_w, span * 0.45)
    if strand == "+":
        body_w = max(span - head_w, scale * 80)
        body = mpatches.Rectangle(
            (left, y - body_h / 2),
            body_w,
            body_h,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=lw,
        )
        tip_x = left + span
        tri = [
            (left + body_w, y - body_h / 2),
            (tip_x, y),
            (left + body_w, y + body_h / 2),
        ]
    else:
        body_w = max(span - head_w, scale * 80)
        body = mpatches.Rectangle(
            (left + head_w, y - body_h / 2),
            body_w,
            body_h,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=lw,
        )
        tip_x = left
        tri = [
            (left + head_w, y - body_h / 2),
            (tip_x, y),
            (left + head_w, y + body_h / 2),
        ]
    head = mpatches.Polygon(
        tri,
        closed=True,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=lw,
    )
    return [body, head]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--delta", type=Path, required=True)
    p.add_argument("--ref-fasta", type=Path, required=True)
    p.add_argument("--qry-fasta", type=Path, required=True)
    p.add_argument("--ref-contig", required=True)
    p.add_argument("--qry-contig", required=True)
    p.add_argument("--gff", type=Path, required=True)
    p.add_argument("--mutations", type=Path, required=True)
    p.add_argument("-o", "--output", type=Path, required=True)
    p.add_argument("--ref-label", default="pre (reference)")
    p.add_argument("--qry-label", default="post (query)")
    p.add_argument(
        "--all-mutation-genes",
        action="store_true",
        help="Label genes from any mutation row on this chromosome (default: unique_to_post only)",
    )
    args = p.parse_args(argv)

    Lr = read_fasta_contig_length(args.ref_fasta, args.ref_contig)
    Lq = read_fasta_contig_length(args.qry_fasta, args.qry_contig)
    segs = parse_show_coords_rows(args.delta, args.ref_contig, args.qry_contig)
    if not segs:
        print(
            f"No show-coords alignments for {args.ref_contig} vs {args.qry_contig}",
            file=sys.stderr,
        )
        return 1

    gff_types = frozenset({"CDS", "gene", "tRNA", "rRNA"})
    genes = parse_gff_genes(args.gff, args.ref_contig, gff_types)
    mut_genes = load_mutation_genes(
        args.mutations, args.ref_contig, only_unique_post=not args.all_mutation_genes
    )

    W = 14.0
    margin = 0.6
    fig, ax = plt.subplots(figsize=(W, 5.2), dpi=150)
    x0 = margin
    plot_w = W - 2 * margin
    sr = plot_w / Lr
    sq = plot_w / Lq

    y_ref = 0.72
    y_qry = 0.28
    bar_h = 0.06

    try:
        cmap = matplotlib.colormaps["tab20"]
    except AttributeError:
        cmap = cm.get_cmap("tab20")
    polys = []
    for i, (r1, r2, q1, q2) in enumerate(segs):
        c = cmap((i % 20) / 19.0)
        color = (*c[:3], 0.45)
        xr1 = x0 + (min(r1, r2) - 1) * sr
        xr2 = x0 + max(r1, r2) * sr
        xq1 = x0 + (min(q1, q2) - 1) * sq
        xq2 = x0 + max(q1, q2) * sq
        poly = mpatches.Polygon(
            [
                (xr1, y_ref + bar_h / 2),
                (xr2, y_ref + bar_h / 2),
                (xq2, y_qry - bar_h / 2),
                (xq1, y_qry - bar_h / 2),
            ],
            closed=True,
            facecolor=color,
            edgecolor="none",
        )
        polys.append(poly)
    for poly in polys:
        ax.add_patch(poly)

    ax.add_patch(
        mpatches.Rectangle(
            (x0, y_ref - bar_h / 2),
            Lr * sr,
            bar_h,
            facecolor="#e8e8e8",
            edgecolor="#333",
            linewidth=0.8,
        )
    )
    ax.add_patch(
        mpatches.Rectangle(
            (x0, y_qry - bar_h / 2),
            Lq * sq,
            bar_h,
            facecolor="#e8e8e8",
            edgecolor="#333",
            linewidth=0.8,
        )
    )

    lbl_patches: list[mpatches.Patch] = []
    label_xy: list[tuple[float, float, str]] = []
    gi = 0
    for gene, hits in sorted(mut_genes.items(), key=lambda kv: min(p for p, _ in kv[1])):
        if gene not in genes:
            continue
        g0, g1, strand = genes[gene]
        center_r = (g0 + g1) / 2.0
        xc_r = x0 + (center_r - 1) * sr
        for patch in gene_arrow_patches(
            g0,
            g1,
            strand,
            y_ref,
            bar_h * 2.8,
            x0,
            sr,
            "#c0392b",
            "#922b21",
            0.6,
        ):
            lbl_patches.append(patch)
        q_center = ref_to_qry(int(round(center_r)), segs)
        if q_center is not None:
            xc_q = x0 + (q_center - 1) * sq
            ax.plot(
                [xc_r, xc_q],
                [y_ref + bar_h * 0.9, y_qry - bar_h * 0.9],
                color="#c0392b",
                lw=0.35,
                alpha=0.5,
                zorder=2,
            )
            tick_h = bar_h * 0.35
            ax.add_patch(
                mpatches.Rectangle(
                    (xc_q - tick_h, y_qry - bar_h / 2 - tick_h),
                    tick_h * 2,
                    tick_h,
                    facecolor="#c0392b",
                    edgecolor="#922b21",
                    linewidth=0.4,
                )
            )
        notes = "; ".join(sorted({n for _, n in hits}))
        short_gene = gene if len(gene) < 22 else gene[:19] + "…"
        txt = f"{short_gene}\n{notes[:70]}{'…' if len(notes) > 70 else ''}"
        stagger = (gi % 5) * 0.034
        label_xy.append((xc_r, y_ref + bar_h * 0.5 + 0.12 + stagger, txt))
        gi += 1

    for patch in lbl_patches:
        patch.set_zorder(4)
        ax.add_patch(patch)

    for xc, yy, txt in label_xy:
        ax.text(
            xc,
            yy,
            txt,
            ha="center",
            va="bottom",
            fontsize=6,
            color="#1a1a1a",
            linespacing=0.95,
            clip_on=False,
        )

    ax.text(
        x0 - 0.08,
        y_ref,
        args.ref_label,
        ha="right",
        va="center",
        fontsize=9,
        fontweight="bold",
    )
    ax.text(
        x0 - 0.08,
        y_qry,
        args.qry_label,
        ha="right",
        va="center",
        fontsize=9,
        fontweight="bold",
    )

    ax.set_xlim(0, W)
    ax.set_ylim(0.08, 1.02)
    ax.axis("off")
    ax.set_title(
        f"Synteny (nucmer): {args.ref_contig} ({Lr:,} bp) vs {args.qry_contig} ({Lq:,} bp)\n"
        "Red arrows: genes with labeled variants"
        + (" (unique to post)" if not args.all_mutation_genes else ""),
        fontsize=10,
        pad=12,
    )

    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
