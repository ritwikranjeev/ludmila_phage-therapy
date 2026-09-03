#!/usr/bin/env python3
"""Genome mutation map for clinical Pre- vs Post-treatment (OCPT04SA01 / OCPT04SA03).

Matches manuscript styling from analyses/Figure5_Phage_resistant_Mutant_Analysis/scripts/plot_strain_domain_mutation_maps.py:
  grey backbone, green mutated segments, two tick rows (pre empty, post variants), same legend layout.
  Title uses the study label "MRSA Clinical Isolate" (not the OCPT strain ID).

Note: This figure uses the MRSA_ONT pre assembly (contig_2). It is not on the same coordinates as
  Figure5 strain_genome_maps/OCPT04SA01_genome_mutation_map.png (S18.OCPT04SA01.Anc); those cannot be merged without liftover.

Inputs (defaults point at MRSA_ONT_comparative outputs):
  - Prokka GenBank for the pre assembly (multi-contig OK; contig_2 is extracted by default)
  - tables/mutations.tsv by default (post-specific SNVs/indels; chrom contig_2). Use
    variants_for_paper.tsv instead for a minimal curated subset (sniff uses column names).

Writes under this analysis folder:
  - figure_inputs/pre_contig2.gbk
  - tables/post_vs_pre_breseq_mutations.tsv
  - figures/MRSA_clinical_isolate_genome_mutation_map*.png|pdf
  - figures/strain_genome_maps_legend.png|pdf (same legend as resistant-mutant maps)
"""

from __future__ import annotations

import csv
import importlib.util
import sys
from collections import defaultdict
from pathlib import Path

# ---- paths -------------------------------------------------------------------
# Supplied by run.sh. This panel reuses the drawing machinery in _strain_maps.py,
# which is located via FIG_SHARED rather than by a relative path.
import os as _os

def _need(v: str) -> Path:
    x = _os.environ.get(v)
    if not x:
        raise SystemExit(f"{v} must be set. Run this figure with ./run.sh")
    return Path(x)

DATA_DIR    = _need("FIG_DATA")
OUT_DIR     = _need("FIG_OUT")
_FIG5_SCRIPT = _need("FIG_SHARED")
OUT_DIR.mkdir(parents=True, exist_ok=True)

DEFAULT_PRE_GBK  = DATA_DIR / "pre.gbk"
DEFAULT_VARIANTS = DATA_DIR / "mutations.tsv"
DEFAULT_TEMPLATE = DATA_DIR / "template_header.tsv"
OUT_INPUT_GBK    = OUT_DIR / "pre_contig2.gbk"
OUT_MUT          = OUT_DIR / "post_vs_pre_breseq_mutations.tsv"
OUT_FIG          = OUT_DIR
LEGEND_BASENAME  = "strain_genome_maps_legend"
CLINICAL_MAP_BASENAME = "MRSA_clinical_isolate_genome_mutation_map"
STRAIN_MAP_FIG_W_IN = 10.0
STRAIN_MAP_FIG_H_IN = 4.35
INTERGENIC_MUT_MARK_HALF_BP = 175.0

def _feature_has_named_gene(feat: dict) -> bool:
    """True if GenBank /gene= is set (exclude locus-only or generic CDS labels)."""
    g = feat.get("gene")
    return bool(g and str(g).strip())


def _load_fig5_module():
    spec = importlib.util.spec_from_file_location("strain_maps", _FIG5_SCRIPT)
    if spec is None or spec.loader is None:
        raise FileNotFoundError(f"Cannot load Figure5 script: {_FIG5_SCRIPT}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["strain_maps"] = mod
    spec.loader.exec_module(mod)
    return mod


def extract_genbank_record(gbk_path: Path, record_contig: str) -> str:
    text = gbk_path.read_text(errors="replace")
    for rec in text.split("\n//\n"):
        rec = rec.strip()
        if not rec:
            continue
        first = rec.split("\n", 1)[0].strip()
        if not first.startswith("LOCUS"):
            continue
        parts = first.split()
        if len(parts) >= 2 and parts[1] == record_contig:
            return rec + "\n//\n"
    raise ValueError(f"No LOCUS record named {record_contig!r} in {gbk_path}")


def breseq_header_from_template(template_tsv: Path) -> list[str]:
    with template_tsv.open() as f:
        return next(csv.reader(f, delimiter="\t"))


def sniff_mutations_pipeline_format(variants_tsv: Path) -> bool:
    """True if MRSA_ONT_comparative-style mutations.tsv (variant_type, unique_to_post)."""
    line = variants_tsv.read_text(errors="replace").splitlines()[0]
    fields = set(line.split("\t"))
    return "variant_type" in fields and "unique_to_post" in fields


def mutations_tsv_to_breseq_rows(
    mutations_tsv: Path,
    chrom: str,
    header: list[str],
    sample_title: str,
) -> list[dict[str, str]]:
    """MRSA_ONT_comparative/tables/mutations.tsv → breseq-compatible rows (post-specific only)."""
    rows_out: list[dict[str, str]] = []
    with mutations_tsv.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if (row.get("chrom") or "").strip() != chrom:
                continue
            up = (row.get("unique_to_post") or "").strip().lower()
            if up not in ("true", "1", "yes"):
                continue
            ref = (row.get("ref") or "").strip()
            alt = (row.get("alt") or "").strip()
            pos = (row.get("pos") or "").strip()
            if not pos:
                continue
            effect = (row.get("effect") or "").strip().lower()
            lr, la = len(ref), len(alt)
            if lr == 1 and la == 1:
                typ = "SNP"
                if "synonymous" in effect:
                    mut_cat = "snp_synonymous"
                else:
                    mut_cat = "snp_nonsynonymous"
            elif lr > la:
                typ = "DEL"
                mut_cat = "deletion"
            elif la > lr:
                typ = "INS"
                mut_cat = "insertion"
            else:
                typ = "SUB"
                mut_cat = "substitution"
            gene = (row.get("gene") or "").strip()
            product = (row.get("product") or "").strip()
            base = {h: "" for h in header}
            base.update(
                {
                    "gene_name": gene,
                    "locus_tag": gene,
                    "gene_product": product,
                    "position": pos,
                    "position_start": pos,
                    "position_end": pos,
                    "mutation_category": mut_cat,
                    "type": typ,
                    "title": sample_title,
                    "ref_seq": ref,
                    "new_seq": alt,
                }
            )
            rows_out.append(base)
    return rows_out


def variants_to_breseq_rows(
    variants_tsv: Path,
    chrom: str,
    header: list[str],
    sample_title: str,
) -> list[dict[str, str]]:
    rows_out: list[dict[str, str]] = []
    with variants_tsv.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if (row.get("chrom") or "").strip() != chrom:
                continue
            ref = (row.get("ref") or "").strip()
            alt = (row.get("alt") or "").strip()
            pos = (row.get("pos") or "").strip()
            if not pos:
                continue
            lr, la = len(ref), len(alt)
            effect = (row.get("effect") or "").strip().lower()
            if lr == 1 and la == 1:
                typ = "SNP"
                if "synonymous" in effect:
                    mut_cat = "snp_synonymous"
                else:
                    mut_cat = "snp_nonsynonymous"
            elif lr > la:
                typ = "DEL"
                mut_cat = "deletion"
            elif la > lr:
                typ = "INS"
                mut_cat = "insertion"
            else:
                typ = "SUB"
                mut_cat = "substitution"
            gene = (row.get("gene") or "").strip()
            lt = (row.get("locus_tag") or "").strip()
            product = (row.get("product") or "").strip()
            base = {h: "" for h in header}
            base.update(
                {
                    "gene_name": gene,
                    "locus_tag": lt or gene,
                    "gene_product": product,
                    "position": pos,
                    "position_start": pos,
                    "position_end": pos,
                    "mutation_category": mut_cat,
                    "type": typ,
                    "title": sample_title,
                    "ref_seq": ref,
                    "new_seq": alt,
                }
            )
            rows_out.append(base)
    return rows_out


def save_legend_like_fig5(sgm) -> None:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    OUT_FIG.mkdir(parents=True, exist_ok=True)
    leg_path = OUT_FIG / f"{LEGEND_BASENAME}.png"
    fig, ax = plt.subplots(figsize=(3.2, 2.2))
    ax.set_axis_off()
    tc = sgm.TYPE_COLORS
    items = [
        *[Patch(facecolor=c, edgecolor="#333333", linewidth=0.5, label=t) for t, c in tc.items()],
        Patch(
            facecolor=sgm.MUTATED_GENE_COLOR,
            edgecolor="#1e6b3a",
            linewidth=0.5,
            label="Mutated gene",
        ),
        Patch(
            facecolor=sgm.TRACK_GREY,
            edgecolor=sgm.TRACK_EDGE,
            linewidth=0.5,
            label="Genome backbone",
        ),
    ]
    ax.legend(
        handles=items,
        loc="center",
        frameon=True,
        fancybox=False,
        edgecolor="#cccccc",
        fontsize=mpl.rcParams["axes.labelsize"],
        title="Genome mutation map",
        title_fontsize=mpl.rcParams["axes.titlesize"],
    )
    fig.tight_layout()
    fig.savefig(leg_path, dpi=300, bbox_inches="tight")
    fig.savefig(OUT_FIG / f"{LEGEND_BASENAME}.pdf", bbox_inches="tight")
    plt.close(fig)


def save_shaped_ticks_legend_like_fig5(sgm) -> None:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    basename = "strain_genome_maps_shaped_ticks_legend"
    leg_path = OUT_FIG / f"{basename}.png"
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    ax.set_axis_off()
    edge = "#2a2a2a"
    items: list = []
    for t, c in sgm.TYPE_COLORS.items():
        if t == "SNP":
            items.append(
                Line2D(
                    [0],
                    [0],
                    marker="|",
                    color=c,
                    markersize=16,
                    markeredgewidth=2.8,
                    markerfacecolor=c,
                    markeredgecolor=c,
                    linestyle="none",
                    label="SNP",
                )
            )
        else:
            mk, ms = sgm.TYPE_MARKERS.get(t, sgm.TYPE_MARKER_DEFAULT)
            items.append(
                Line2D(
                    [0],
                    [0],
                    marker=mk,
                    color="w",
                    markerfacecolor=c,
                    markeredgecolor=edge,
                    markeredgewidth=0.55,
                    markersize=ms,
                    linestyle="none",
                    label=t,
                )
            )
    items.append(
        Patch(
            facecolor=sgm.MUTATED_GENE_COLOR,
            edgecolor="#1e6b3a",
            linewidth=0.5,
            label="Mutated gene",
        )
    )
    items.append(
        Patch(
            facecolor=sgm.TRACK_GREY,
            edgecolor=sgm.TRACK_EDGE,
            linewidth=0.5,
            label="Genome backbone",
        )
    )
    ax.legend(
        handles=items,
        loc="center",
        frameon=True,
        fancybox=False,
        edgecolor="#cccccc",
        fontsize=mpl.rcParams["axes.labelsize"],
        title="Genome map (shaped ticks)",
        title_fontsize=mpl.rcParams["axes.titlesize"],
    )
    fig.tight_layout()
    fig.savefig(leg_path, dpi=300, bbox_inches="tight")
    fig.savefig(OUT_FIG / f"{basename}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_clinical_genome_map(
    sgm,
    ref_gbk: Path,
    mut_tsv: Path,
    *,
    shaped_mutation_markers: bool,
) -> Path:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    sgm.manuscript_rc()

    genome_len, cds_list = sgm.parse_cds_features(ref_gbk)
    if genome_len <= 0:
        raise RuntimeError(f"No genome length from {ref_gbk}")

    locus_to_sseqid: dict[str, str] = {}

    def load_mut_local(path: Path) -> list[dict]:
        if not path.is_file():
            return []
        with path.open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    rep_rows = {
        1: [],  # Pre-treatment: no variants vs self on this coordinate system
        2: load_mut_local(mut_tsv),  # Post-treatment vs pre
    }
    mutated = sgm.mutated_gene_names(rep_rows)
    mutated_feats = [f for f in cds_list if sgm.feature_is_mutated(f, mutated)]
    mutated_sorted = sorted(mutated_feats, key=lambda f: (f["start"], f["end"]))

    label_xm: list[float] = []
    label_texts: list[str] = []
    for feat in mutated_sorted:
        if sgm.mutated_feat_has_diamond_hit(feat, locus_to_sseqid):
            continue
        if not _feature_has_named_gene(feat):
            continue
        label_xm.append((feat["start"] + feat["end"]) / 2.0)
        label_texts.append(sgm.label_mutated_cds(feat))

    layout = sgm.layout_mutated_gene_labels(
        label_xm,
        label_texts,
        float(genome_len),
        STRAIN_MAP_FIG_W_IN,
        sgm.GENE_LABEL_FONTSIZE_PT,
    )
    x_placed = [layout[i][0] for i in range(len(layout))]
    row_placed = [layout[i][1] for i in range(len(layout))]
    lane_h = sgm.assign_leader_horizontal_lanes(
        label_xm, x_placed, row_placed, float(genome_len)
    )
    text_stagger = sgm.assign_label_text_stagger(
        layout,
        lane_h,
        sgm.LABEL_TEXT_STAGGER_TIERS,
        sgm.LABEL_TEXT_STAGGER_DY,
    )

    y_track = 0.0
    track_half = 0.30
    # Same vertical ordering as Figure5 replicates: rep 1 higher on the axis
    y_reps = {1: 2.0, 2: 1.0}
    reps_to_plot = [1, 2]
    tick_half = 0.24 if shaped_mutation_markers else 0.2

    y_bar_bottom = y_track - track_half
    leader_stub_drop = 0.28
    label_row_step = 0.52
    text_below_arm = 0.26
    leader_lane_dy = 0.062

    fig, ax = plt.subplots(figsize=(STRAIN_MAP_FIG_W_IN, STRAIN_MAP_FIG_H_IN))

    ax.barh(
        y_track,
        genome_len,
        left=0,
        height=2 * track_half,
        color=sgm.TRACK_GREY,
        edgecolor=sgm.TRACK_EDGE,
        linewidth=0.4,
        zorder=2,
    )

    for feat in mutated_sorted:
        s, e = feat["start"], feat["end"]
        width = e - s + 1
        ax.barh(
            y_track,
            width,
            left=s,
            height=2 * track_half,
            color=sgm.MUTATED_GENE_COLOR,
            edgecolor="#1e6b3a",
            linewidth=0.35,
            zorder=3,
        )

    for p in sgm.intergenic_mutation_positions(rep_rows, cds_list, float(genome_len)):
        x0 = max(0.0, float(p) - INTERGENIC_MUT_MARK_HALF_BP)
        x1 = min(float(genome_len), float(p) + INTERGENIC_MUT_MARK_HALF_BP)
        w = x1 - x0
        if w < 1.0:
            continue
        ax.barh(
            y_track,
            w,
            left=x0,
            height=2 * track_half,
            color=sgm.MUTATED_GENE_COLOR,
            edgecolor="#1e6b3a",
            linewidth=0.35,
            zorder=3,
        )

    for i in range(len(label_xm)):
        xm_gene = label_xm[i]
        lab = label_texts[i]
        x_label, row = layout[i]
        lh = lane_h[i]
        y_arm = y_bar_bottom - leader_stub_drop - row * label_row_step - lh * leader_lane_dy
        y_text = y_arm - text_below_arm - text_stagger[i]
        ax.plot(
            [xm_gene, xm_gene, x_label, x_label],
            [y_bar_bottom, y_arm, y_arm, y_text],
            linestyle="-",
            color="#6e6e6e",
            linewidth=0.65,
            alpha=0.88,
            clip_on=False,
            zorder=2.5,
            solid_capstyle="round",
        )
        ax.text(
            x_label,
            y_text,
            lab,
            rotation=0,
            va="top",
            ha="center",
            fontsize=sgm.GENE_LABEL_FONTSIZE_PT,
            color="#222222",
            clip_on=False,
            zorder=4,
        )

    y_track_top = y_track + track_half

    for rep in reps_to_plot:
        y = y_reps[rep]
        ax.axhline(y, xmin=0, xmax=1, color="#e8e8e8", linewidth=0.75, zorder=1)
        rows = rep_rows.get(rep, [])

        if shaped_mutation_markers:
            events: list[tuple[int, str]] = []
            for r in rows:
                p = sgm.as_int_pos(r)
                if p is None or p < 1 or p > genome_len:
                    continue
                t = (r.get("type") or "SNP").strip().upper()
                events.append((p, t))
            for x_draw, t, _bp in sgm.cluster_mutation_x_offsets(events, float(genome_len)):
                color = sgm.TYPE_COLORS.get(t, "#4d4d4d")
                ax.plot(
                    [x_draw, x_draw],
                    [y_track_top, y - tick_half * 0.9],
                    linestyle=":",
                    color=color,
                    linewidth=1.05,
                    alpha=0.52,
                    solid_capstyle="round",
                    zorder=3,
                )
                if t == "SNP":
                    ax.vlines(
                        x_draw,
                        y - tick_half,
                        y + tick_half,
                        color=color,
                        linewidth=2.0,
                        alpha=0.96,
                        zorder=4,
                    )
                else:
                    marker, ms = sgm.TYPE_MARKERS.get(t, sgm.TYPE_MARKER_DEFAULT)
                    ax.plot(
                        x_draw,
                        y,
                        marker=marker,
                        ms=ms,
                        color=color,
                        markeredgecolor="#2a2a2a",
                        markeredgewidth=0.55,
                        linestyle="none",
                        zorder=5,
                        clip_on=False,
                    )
        else:
            pos_by_type: dict[str, list[int]] = defaultdict(list)
            for r in rows:
                p = sgm.as_int_pos(r)
                if p is None or p < 1 or p > genome_len:
                    continue
                t = (r.get("type") or "SNP").strip().upper()
                pos_by_type[t].append(p)
            for t, positions in pos_by_type.items():
                color = sgm.TYPE_COLORS.get(t, "#4d4d4d")
                y_tick_bottom = y - tick_half
                for p in positions:
                    ax.plot(
                        [p, p],
                        [y_track_top, y_tick_bottom],
                        linestyle=":",
                        color=color,
                        linewidth=1.05,
                        alpha=0.52,
                        solid_capstyle="round",
                        zorder=3,
                    )
                ax.vlines(
                    positions,
                    y - tick_half,
                    y + tick_half,
                    color=color,
                    linewidth=2.0,
                    alpha=0.96,
                    zorder=4,
                )

    if label_xm:
        lowest_text_y = y_bar_bottom
        for i in range(len(label_xm)):
            row = layout[i][1]
            lh = lane_h[i]
            y_arm = y_bar_bottom - leader_stub_drop - row * label_row_step - lh * leader_lane_dy
            lowest_text_y = min(lowest_text_y, y_arm - text_below_arm - text_stagger[i])
        label_depth = (y_bar_bottom - lowest_text_y) + 0.34
    else:
        label_depth = 0.35

    if shaped_mutation_markers:
        x_pad = max(5000.0, float(genome_len) / 2200.0)
        ax.set_xlim(0, genome_len * 1.002 + 0.4 * x_pad)
    else:
        ax.set_xlim(0, genome_len * 1.002)

    y_top = 2.62
    ax.set_ylim(y_track - track_half - label_depth, y_top)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _p: f"{x / 1e6:.2f}"))
    ax.set_yticks([])
    ax.tick_params(axis="both", labelsize=mpl.rcParams["xtick.labelsize"])
    ax.tick_params(axis="y", length=0, labelleft=False)
    ax.set_title(
        "MRSA Clinical Isolate",
        fontweight="bold",
        pad=14,
    )
    ax.grid(axis="x", color="#f0f0f0", linewidth=0.5)

    OUT_FIG.mkdir(parents=True, exist_ok=True)
    suffix = "_shaped_ticks" if shaped_mutation_markers else ""
    out = OUT_FIG / f"{CLINICAL_MAP_BASENAME}{suffix}.png"
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.20, left=0.10, top=0.92)
    fig.savefig(out, dpi=300, bbox_inches=None, pad_inches=0)
    fig.savefig(out.with_suffix(".pdf"), bbox_inches=None, pad_inches=0)
    plt.close(fig)
    return out


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pre-gbk", type=Path, default=DEFAULT_PRE_GBK)
    ap.add_argument(
        "--variants",
        type=Path,
        default=DEFAULT_VARIANTS,
        help="MRSA_ONT mutations.tsv (default) or variants_for_paper.tsv",
    )
    ap.add_argument("--chrom", default="contig_2")
    ap.add_argument(
        "--template-mutations",
        type=Path,
        default=DEFAULT_TEMPLATE,
        help="Existing breseq-style TSV whose header row is copied for output compatibility",
    )
    args = ap.parse_args()

    if not args.pre_gbk.is_file():
        raise SystemExit(f"Missing pre GenBank: {args.pre_gbk}")
    if not args.variants.is_file():
        raise SystemExit(f"Missing variants TSV: {args.variants}")
    if not args.template_mutations.is_file():
        raise SystemExit(f"Missing template mutations TSV: {args.template_mutations}")

    sgm = _load_fig5_module()

    OUT_INPUT_GBK.parent.mkdir(parents=True, exist_ok=True)
    rec = extract_genbank_record(args.pre_gbk, args.chrom)
    OUT_INPUT_GBK.write_text(rec)
    print(f"Wrote {OUT_INPUT_GBK}")

    header = breseq_header_from_template(args.template_mutations)
    if sniff_mutations_pipeline_format(args.variants):
        rows = mutations_tsv_to_breseq_rows(
            args.variants,
            args.chrom,
            header,
            "OCPT04SA03_vs_OCPT04SA01",
        )
    else:
        rows = variants_to_breseq_rows(
            args.variants,
            args.chrom,
            header,
            "OCPT04SA03_vs_OCPT04SA01",
        )
    OUT_MUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT_MUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)
    print(f"Wrote {OUT_MUT} ({len(rows)} variants)")

    save_legend_like_fig5(sgm)
    print(f"Wrote {OUT_FIG / (LEGEND_BASENAME + '.png')}")
    save_shaped_ticks_legend_like_fig5(sgm)
    print(f"Wrote {OUT_FIG / 'strain_genome_maps_shaped_ticks_legend.png'}")

    p1 = plot_clinical_genome_map(sgm, OUT_INPUT_GBK, OUT_MUT, shaped_mutation_markers=False)
    print(f"Wrote {p1}")
    p2 = plot_clinical_genome_map(sgm, OUT_INPUT_GBK, OUT_MUT, shaped_mutation_markers=True)
    print(f"Wrote {p2}")


if __name__ == "__main__":
    main()
