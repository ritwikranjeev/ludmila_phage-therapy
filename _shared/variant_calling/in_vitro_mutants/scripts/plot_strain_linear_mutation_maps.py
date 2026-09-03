#!/usr/bin/env python3
"""Per-strain linear genome mutation maps with replicate tracks.

Creates one plot per host strain:
- X-axis: linear genome coordinate
- Five horizontal replicate tracks above baseline (rep 1..5)
- Mutation marks per replicate
- "Important genes" labeled on baseline
"""

from __future__ import annotations

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
SAMPLES_TSV = ROOT / "metadata" / "samples.tsv"
MUT_DIR = ROOT / "results" / "tables" / "mutations"
REF_DIR = ROOT / "results" / "references"
OUT_DIR = ROOT / "results" / "figures" / "strain_linear_maps"

TYPE_COLORS = {
    "SNP": "#2e86ab",
    "SUB": "#8e44ad",
    "DEL": "#f18f01",
    "INS": "#bc4b51",
}


def parse_gbk(gbk_path: Path) -> tuple[int, dict[str, tuple[int, int]]]:
    """Return genome length and gene -> (start, end)."""
    genome_len = 0
    genes: dict[str, tuple[int, int]] = {}

    current_span = None
    current_gene = None
    current_locus = None
    in_features = False

    with gbk_path.open() as handle:
        for line in handle:
            if line.startswith("LOCUS"):
                m = re.search(r"\b(\d+)\s+bp\b", line)
                if m:
                    genome_len = int(m.group(1))
            if line.startswith("FEATURES"):
                in_features = True
                continue
            if not in_features:
                continue
            if line.startswith("ORIGIN"):
                break

            if line.startswith("     CDS"):
                # flush previous CDS
                if current_span and (current_gene or current_locus):
                    genes[current_gene or current_locus] = current_span
                current_gene = None
                current_locus = None
                nums = [int(x) for x in re.findall(r"\d+", line)]
                if nums:
                    current_span = (min(nums), max(nums))
                else:
                    current_span = None
                continue

            if current_span is None:
                continue

            if '/gene="' in line:
                m = re.search(r'/gene="([^"]+)"', line)
                if m:
                    current_gene = m.group(1)
            elif '/locus_tag="' in line:
                m = re.search(r'/locus_tag="([^"]+)"', line)
                if m:
                    current_locus = m.group(1)

    if current_span and (current_gene or current_locus):
        genes[current_gene or current_locus] = current_span

    return genome_len, genes


def load_samples() -> dict[str, dict]:
    out = {}
    with SAMPLES_TSV.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            out[row["sample_id"]] = row
    return out


def load_mutations(sample_id: str) -> list[dict]:
    path = MUT_DIR / f"{sample_id}.tsv"
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def choose_important_genes(rep_rows: dict[int, list[dict]], top_n: int = 8) -> list[str]:
    by_gene_reps = defaultdict(set)
    by_gene_events = Counter()
    for rep, rows in rep_rows.items():
        for r in rows:
            g = (r.get("gene_name") or "").strip()
            if not g:
                continue
            # split composite names like arlS/odhA
            for part in g.split("/"):
                gene = part.strip()
                if gene:
                    by_gene_reps[gene].add(rep)
                    by_gene_events[gene] += 1
    ranked = sorted(
        by_gene_reps.keys(),
        key=lambda g: (-len(by_gene_reps[g]), -by_gene_events[g], g),
    )
    return ranked[:top_n]


def merge_intervals(intervals: list[tuple[int, int]], gap: int = 300) -> list[tuple[int, int]]:
    """Merge close/overlapping genomic intervals for cleaner genome band display."""
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        ms, me = merged[-1]
        if s <= me + gap:
            merged[-1] = (ms, max(me, e))
        else:
            merged.append((s, e))
    return merged


def plot_strain(host_id: str, rep_to_sample: dict[int, str], ref_gbk: Path) -> Path:
    genome_len, gene_coords = parse_gbk(ref_gbk)

    rep_rows: dict[int, list[dict]] = {rep: load_mutations(sid) for rep, sid in rep_to_sample.items()}
    important = choose_important_genes(rep_rows, top_n=10)
    mutated_genes = set()
    for rows in rep_rows.values():
        for r in rows:
            g = (r.get("gene_name") or "").strip()
            if not g:
                continue
            for part in g.split("/"):
                part = part.strip()
                if part:
                    mutated_genes.add(part)
    mutated_spans = [gene_coords[g] for g in mutated_genes if g in gene_coords]
    mutated_band_spans = merge_intervals(mutated_spans, gap=250)

    fig, ax = plt.subplots(figsize=(14, 5))
    y_base = 0
    y_tracks = {1: 5, 2: 4, 3: 3, 4: 2, 5: 1}

    # baseline genome + fat band where mutated genes are located
    ax.hlines(y=y_base, xmin=0, xmax=genome_len, color="#222222", linewidth=2, zorder=1)
    for s, e in mutated_band_spans:
        ax.hlines(y=y_base, xmin=s, xmax=e, color="#2ecc71", linewidth=8, alpha=0.65, zorder=2)

    # important gene spans + labels
    label_y = -0.55
    for idx, g in enumerate(important):
        span = gene_coords.get(g)
        if not span:
            continue
        s, e = span
        ax.hlines(y=-0.08, xmin=s, xmax=e, color="#7f8c8d", linewidth=3, alpha=0.8)
        xmid = (s + e) / 2
        ax.text(
            xmid,
            label_y - 0.35 * (idx % 2),
            g,
            fontsize=8,
            rotation=45,
            ha="left",
            va="top",
        )

    # replicate tracks
    for rep in range(1, 6):
        y = y_tracks[rep]
        sid = rep_to_sample.get(rep, f"rep{rep}")
        rows = rep_rows.get(rep, [])
        ax.hlines(y=y, xmin=0, xmax=genome_len, color="#d9d9d9", linewidth=1)

        if rows:
            pos_by_type = defaultdict(list)
            for r in rows:
                try:
                    p = int(float(r.get("position", "nan")))
                except Exception:
                    continue
                t = (r.get("type") or "SNP").strip().upper()
                pos_by_type[t].append(p)
            for t, positions in pos_by_type.items():
                color = TYPE_COLORS.get(t, "#4d4d4d")
                ax.vlines(positions, y - 0.35, y + 0.35, color=color, linewidth=0.9, alpha=0.9)
            ax.text(-genome_len * 0.01, y, f"R{rep}: {sid} (n={len(rows)})", ha="right", va="center", fontsize=8)
        else:
            ax.text(-genome_len * 0.01, y, f"R{rep}: {sid} (missing)", ha="right", va="center", fontsize=8, color="#9e9e9e")

    ax.set_xlim(0, genome_len)
    ax.set_ylim(-1.4, 5.8)
    ax.set_xlabel("Genome position (bp)")
    ax.set_yticks([])
    ax.set_title(f"{host_id}: replicate mutation tracks on linear reference genome")

    handles = [
        plt.Line2D([0], [0], color=c, lw=2, label=t)
        for t, c in TYPE_COLORS.items()
    ]
    handles.append(plt.Line2D([0], [0], color="#2ecc71", lw=6, label="mutated-gene regions"))
    ax.legend(handles=handles, title="Mutation type", loc="upper right", frameon=False)
    ax.grid(axis="x", color="#efefef", linewidth=0.6)
    fig.tight_layout()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / f"{host_id}_linear_mutation_map.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    return out


def main() -> None:
    samples = load_samples()

    # group resistant mutants by host_id and replicate
    groups: dict[str, dict[int, str]] = defaultdict(dict)
    ref_by_host: dict[str, str] = {}
    for sid, row in samples.items():
        if row["sample_type"] == "host_reference":
            ref_by_host[row["host_id"]] = row["sample_id"]
        elif row["sample_type"] == "resistant_mutant":
            rep_raw = row.get("replicate", "").strip()
            if rep_raw.isdigit():
                groups[row["host_id"]][int(rep_raw)] = sid

    outputs = []
    for host_id, rep_map in sorted(groups.items()):
        ref_id = ref_by_host.get(host_id)
        if not ref_id:
            continue
        gbk = REF_DIR / ref_id / f"{ref_id}.gbk"
        if not gbk.exists():
            continue
        outputs.append(plot_strain(host_id, rep_map, gbk))

    print("Generated:")
    for p in outputs:
        print(f"- {p}")


if __name__ == "__main__":
    main()
