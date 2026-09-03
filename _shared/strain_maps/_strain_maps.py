#!/usr/bin/env python3
"""Full-genome mutation maps per strain (manuscript-style).

- Single linear genome on x-axis (bp)
- One continuous grey genome backbone; mutated CDS, rRNA, tRNA, etc. drawn in green on top; intergenic mutations as short green segments on the backbone
- Labels only under mutated (green) segments; every green segment gets a label (GBK-only; no label if DIAMOND hit). OCPT04SA01 main maps: skip text labels for features without /gene= (locus-only / generic CDS labels suppressed)
- Labels are horizontal, packed into rows with solid L-shaped leaders (GBK-only; no label if DIAMOND hit)
- Rep. 1–5 rows with short mutation ticks (hypermutator replicates excluded from main maps); OCPT04SA01 adds rep. 6 when results/tables/mutations/MRSA_clinical_lifted_S18_chr1.tsv exists (post-treatment clinical SNVs/indels lifted onto S18 chr1)
- Main maps: “Hypermutator” badge at top of axes (strains with a hypermutator panel); OCPT04SA01 rep. 4 shows “Excluded” in-track (grey box); replicate tag font sizes scale from y-axis rep label size (smaller)
- Separate one-row maps in strain_genome_maps_hypermutators/; mutation ticks use horizontal cluster offsets (same as shaped-tick maps) so dense variants are not drawn on one vertical line; OCPT04SA01 labels = post-treatment (main-map) genes plus phage host-factor CDS (see feat_matches_phage_host_factor); PedCF hypermutator x-axis zooms to point-variant span (megabase structural rows excluded from zoom extent)
- Standalone legend: strain_hypermutator_maps_legend.png (hypermutator folder only)
- Optional duplicate maps with shaped markers per variant type in strain_genome_maps_shaped_ticks/
- Legend in a separate figure (and a second legend for shaped ticks)
"""

from __future__ import annotations

import csv
import io
import re
import ssl
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter

# ---- paths -------------------------------------------------------------------
# Supplied by run.sh. Nothing here knows where it sits in the tree.
import os as _os

def _need(v: str) -> Path:
    x = _os.environ.get(v)
    if not x:
        raise SystemExit(f"{v} must be set. Run this figure with ./run.sh")
    return Path(x)

DATA_DIR = _need("FIG_DATA")
OUT_ROOT = _need("FIG_OUT")

SAMPLES_TSV = DATA_DIR / "samples.tsv"
MUT_DIR     = DATA_DIR / "mutations"
REF_DIR     = DATA_DIR / "references"
DIAMOND_DIR = DATA_DIR / "diamond"

OUT_DIR              = OUT_ROOT
OUT_DIR_SHAPED       = OUT_ROOT / "shaped_ticks"
OUT_DIR_HYPERMUTATOR = OUT_ROOT / "hypermutators"

OCPT04SA01_CLINICAL_MUT_TSV = MUT_DIR / "MRSA_clinical_lifted_S18_chr1.tsv"
OCPT04SA01_CLINICAL_REP = 6
LEGEND_BASENAME = "strain_genome_maps_legend"
LEGEND_SHAPED_BASENAME = "strain_genome_maps_shaped_ticks_legend"
LEGEND_HYPERMUTATOR_BASENAME = "strain_hypermutator_maps_legend"
UNIPROT_CACHE = DIAMOND_DIR / "uniprot_label_cache.tsv"

# Only render the strains named here (comma-separated). Unset = all.
ONLY_STRAINS = [x for x in _os.environ.get("FIG_STRAINS", "").split(",") if x]

# Which map types to render (comma-separated): main, shaped, hypermutator.
# Unset = all three. One strain feeds several display items — the OCPT04SA01
# main map is Figure 5a while its hypermutator map is Supplementary 5d — so the
# caller selects the map rather than getting every map the strain can produce.
_MAP_KINDS = ("main", "shaped", "hypermutator")
ONLY_MAPS = [x.strip() for x in _os.environ.get("FIG_MAPS", "").split(",") if x.strip()]
for _m in ONLY_MAPS:
    if _m not in _MAP_KINDS:
        raise SystemExit(f"FIG_MAPS: unknown map type {_m!r}; choose from {', '.join(_MAP_KINDS)}")
if not ONLY_MAPS:
    ONLY_MAPS = list(_MAP_KINDS)


def want_map(kind: str) -> bool:
    return kind in ONLY_MAPS

# Network access to UniProt is opt-in; the committed cache is used otherwise.
ALLOW_UNIPROT_FETCH = _os.environ.get("FIG_ALLOW_NETWORK", "0") == "1"

UNIPROT_CHUNK = 35
UNIPROT_SLEEP_S = 0.35


def _uniprot_ssl_context() -> ssl.SSLContext:
    try:
        import certifi

        return ssl.create_default_context(cafile=certifi.where())
    except Exception:
        return ssl.create_default_context()


def urlopen_uniprot(req: urllib.request.Request, timeout: float = 120):
    if not ALLOW_UNIPROT_FETCH:
        raise RuntimeError(
            "UniProt lookup attempted but FIG_ALLOW_NETWORK is not set. "
            "The committed cache should cover all accessions; set FIG_ALLOW_NETWORK=1 "
            "only if you intend to refresh it."
        )
    """HTTPS GET with certs; on macOS/Framework Python cert issues, retry without verify."""
    try:
        return urllib.request.urlopen(req, timeout=timeout, context=_uniprot_ssl_context())
    except urllib.error.URLError as e:
        msg = str(e).lower()
        if "certificate verify failed" in msg:
            ctx = ssl._create_unverified_context()
            return urllib.request.urlopen(req, timeout=timeout, context=ctx)
        raise

HYPERMUTATOR_SAMPLES = {"S15.OCPT04SA01.3.6", "S26.PedCF.2.6"}

# No mutation ticks / empty rep row on multi-replicate main maps (supplementary single-rep plot instead)
HYPERMUTATOR_EXCLUDE_FROM_MAIN_MAP = frozenset(
    {"S15.OCPT04SA01.3.6", "S26.PedCF.2.6"}
)

# Host → replicate for hypermutator-only supplementary maps
HYPERMUTATOR_REP_BY_HOST: dict[str, int] = {
    "OCPT04SA01": 3,
    "PedCF": 2,
}

# De-emphasize this replicate row (still shows true positions; visual only)
GREY_OUT_REP4_HOSTS = {"OCPT04SA01"}

GENE_LABEL_FONTSIZE_PT = 10.5
# Extra vertical drop for label text (tiers cycle left→right per leader row/lane)
LABEL_TEXT_STAGGER_TIERS = 6
LABEL_TEXT_STAGGER_DY = 0.14
# x_center past this fraction of genome length: pack labels by trying leftward offsets first
RIGHT_PACK_GENOME_FRAC = 0.54

TYPE_COLORS = {
    "SNP": "#2e86ab",
    "SUB": "#8e44ad",
    "DEL": "#f18f01",
    "INS": "#bc4b51",
}

# De-emphasized row: same hue families, lower chroma (positions and classes unchanged)
TYPE_COLORS_MUTED = {
    "SNP": "#6f92a8",
    "SUB": "#8f7aa3",
    "DEL": "#c49a62",
    "INS": "#a77b7f",
}

# Matplotlib marker name, markersize (points) for shaped-tick variant (SNP uses bar ticks, not these)
TYPE_MARKERS: dict[str, tuple[str, float]] = {
    "SUB": ("s", 7.2),
    "DEL": ("v", 8.0),
    "INS": ("D", 7.5),
}
TYPE_MARKER_DEFAULT: tuple[str, float] = ("^", 6.8)

MUTATED_GENE_COLOR = "#27ae60"
TRACK_GREY = "#d8d8d8"
TRACK_EDGE = "#b0b0b0"

# Figure size (inches) for strain maps
STRAIN_MAP_FIG_W_IN = 10.0
STRAIN_MAP_FIG_H_IN = 4.0
STRAIN_MAP_FIG_H_HYPERMUTATOR_IN = 2.75
# OCPT04SA01 hypermutator: extra label density (smaller text, wider fig, looser rows)
OCPT_HM_LABEL_FONT_SCALE = 0.82
OCPT_HM_LAYOUT_GAP_MULT = 1.55
OCPT_HM_LAYOUT_STEP_MULT = 0.88
OCPT_HM_LABEL_ROW_STEP = 0.58
OCPT_HM_FIG_W_MULT = 1.14
OCPT_HM_FIG_H_MULT = 1.28

# PedCF hypermutator zoom: ensure x-axis extends to at least this position (Mb) on the right
PEDCF_HYPERMUTATOR_XLIM_RIGHT_MIN_MB = 0.6

# Half-width (bp) for intergenic mutation highlights (no annotated feature interval)
INTERGENIC_MUT_MARK_HALF_BP = 175.0


def manuscript_rc() -> None:
    mpl.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "font.family": "serif",
            "font.serif": [
                "Times New Roman",
                "Times",
                "Nimbus Roman",
                "DejaVu Serif",
                "serif",
            ],
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 12,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


# In-track “Hypermutator” / “Excluded” tags (smaller than y-axis rep numbers)
ROW_REPLICATE_TAG_FONT_SCALE = 0.52


def _replicate_row_tag_fontsize() -> float:
    return float(mpl.rcParams["ytick.labelsize"]) * ROW_REPLICATE_TAG_FONT_SCALE


def draw_hypermutator_axes_badge(ax) -> None:
    """Top-of-panel label for OCPT04SA01 / PedCF maps (main and hypermutator-only)."""
    ax.text(
        0.99,
        0.995,
        "Hypermutator",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=float(mpl.rcParams["axes.titlesize"]) * 0.92,
        fontweight="bold",
        color="#1a1a1a",
        bbox=dict(
            boxstyle="round,pad=0.35",
            facecolor="#fff8e8",
            edgecolor="#c9a227",
            linewidth=0.85,
        ),
        zorder=15,
        clip_on=False,
    )


def draw_excluded_replicate_tag(ax, genome_len: float, y: float) -> None:
    ax.text(
        float(genome_len) * 0.992,
        y,
        "Excluded",
        ha="right",
        va="center",
        fontsize=_replicate_row_tag_fontsize(),
        color="#3d3d3d",
        bbox=dict(
            boxstyle="round,pad=0.32",
            facecolor="#e4e4e4",
            edgecolor="#9c9c9c",
            linewidth=0.75,
        ),
        zorder=6,
        clip_on=False,
    )


def _parse_product(chunk: str) -> str | None:
    if '/product="' not in chunk:
        return None
    lines = chunk.split("\n")
    for i, line in enumerate(lines):
        if '/product="' not in line:
            continue
        after = line.split('/product="', 1)[1]
        parts: list[str] = []
        if '"' in after:
            parts.append(after.split('"', 1)[0])
            return " ".join(" ".join(parts).split()).strip() or None
        parts.append(after.strip())
        for cont in lines[i + 1 :]:
            t = cont.strip()
            if not t:
                continue
            if t.endswith('"'):
                parts.append(t[:-1].strip())
                return " ".join(" ".join(parts).split()).strip() or None
            parts.append(t)
        return None
    return None


_GBK_FEATURE_LINE = re.compile(
    r"\n     (CDS|rRNA|tRNA|tmRNA|misc_RNA)\s+", re.MULTILINE
)


def parse_cds_features(gbk_path: Path) -> tuple[int, list[dict]]:
    """Parse CDS and non-coding RNA features from GenBank (shared coordinates/gene fields)."""
    text = gbk_path.read_text(errors="replace")
    glen = 0
    m = re.search(r"LOCUS\s+[^\n]+\s+(\d+)\s+bp\b", text)
    if m:
        glen = int(m.group(1))

    feats: list[dict] = []
    for mfeat in _GBK_FEATURE_LINE.finditer(text):
        start_ch = mfeat.end()
        next_m = _GBK_FEATURE_LINE.search(text, mfeat.end())
        end_ch = next_m.start() if next_m else len(text)
        chunk = text[start_ch:end_ch]
        first = chunk.split("\n", 1)[0]
        nums = [int(x) for x in re.findall(r"\d+", first)]
        if len(nums) < 2:
            continue
        start, end = min(nums), max(nums)
        gene_m = re.search(r'/gene="([^"]+)"', chunk)
        locus_m = re.search(r'/locus_tag="([^"]+)"', chunk)
        product = _parse_product(chunk)
        feats.append(
            {
                "start": start,
                "end": end,
                "gene": gene_m.group(1).strip() if gene_m else None,
                "locus_tag": locus_m.group(1).strip() if locus_m else None,
                "product": product,
            }
        )
    if glen <= 0 and feats:
        glen = max(f["end"] for f in feats)
    return glen, feats


def display_label(feat: dict) -> str | None:
    g = feat.get("gene")
    if g and g.strip():
        return g.strip()
    prod = feat.get("product") or ""
    prod = prod.strip()
    if not prod:
        return None
    low = prod.lower()
    if "hypothetical" in low and "protein" in low:
        return None
    if len(prod) > 42:
        return prod[:39].rstrip() + "…"
    return prod


def label_mutated_cds(feat: dict) -> str:
    """Always produce a label for a mutated (green) CDS."""
    s = display_label(feat)
    if s:
        return s
    lt = feat.get("locus_tag")
    if lt and lt.strip():
        return lt.strip()
    prod = (feat.get("product") or "").strip()
    if prod:
        return (prod[:40] + "…") if len(prod) > 40 else prod
    return f"CDS {feat['start']}–{feat['end']}"


def feature_has_named_gene(feat: dict) -> bool:
    """True if GenBank /gene= is set (skip locus-only / generic CDS labels)."""
    g = feat.get("gene")
    return bool(g and str(g).strip())


def load_diamond_qseqid_map(tsv_path: Path) -> dict[str, str]:
    """Map query locus (qseqid) -> subject sseqid from DIAMOND best-hit TSV."""
    if not tsv_path.is_file():
        return {}
    out: dict[str, str] = {}
    with tsv_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            q = (row.get("qseqid") or "").strip()
            s = (row.get("sseqid") or "").strip()
            if q and s:
                out.setdefault(q, s)
    return out


UniProtMeta = tuple[str, str, str]  # gene_primary, gene_names, protein_name


def _uniprot_row_nonempty(meta: UniProtMeta) -> bool:
    return any((x or "").strip() for x in meta)


def load_uniprot_cache() -> dict[str, UniProtMeta]:
    """accession -> (gene_primary, gene_names cell, protein_name). Supports legacy 3-column cache."""
    if not UNIPROT_CACHE.is_file():
        return {}
    d: dict[str, UniProtMeta] = {}
    with UNIPROT_CACHE.open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows:
        return {}
    header = [h.strip() for h in rows[0]]
    hmap = {h: i for i, h in enumerate(header)}
    i_entry = hmap.get("Entry", 0)
    i_gp = hmap.get("Gene Names (primary)")
    i_gn = hmap.get("Gene Names")
    i_pn = hmap.get("Protein names")
    # Legacy header: Entry, Gene Names, Protein names
    if i_gn is None and len(header) >= 3:
        i_gn, i_pn = 1, 2
    for row in rows[1:]:
        if not row or i_entry >= len(row):
            continue
        acc = row[i_entry].strip()
        if not acc:
            continue
        gp = row[i_gp].strip() if i_gp is not None and i_gp < len(row) else ""
        gn = row[i_gn].strip() if i_gn is not None and i_gn < len(row) else ""
        pn = row[i_pn].strip() if i_pn is not None and i_pn < len(row) else ""
        d[acc] = (gp, gn, pn)
    return d


def write_uniprot_cache(cache: dict[str, UniProtMeta]) -> None:
    UNIPROT_CACHE.parent.mkdir(parents=True, exist_ok=True)
    with UNIPROT_CACHE.open("w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t")
        w.writerow(["Entry", "Gene Names (primary)", "Gene Names", "Protein names"])
        for acc in sorted(cache.keys()):
            gp, gn, pn = cache[acc]
            w.writerow([acc, gp, gn, pn])


def fetch_uniprot_chunk(accessions: list[str]) -> dict[str, UniProtMeta]:
    """GET UniProt TSV stream for a batch of accessions."""
    if not accessions:
        return {}
    inner = ") OR (".join(f"accession:{a}" for a in accessions)
    query = f"({inner})"
    params = urllib.parse.urlencode(
        {
            "format": "tsv",
            "fields": "accession,gene_primary,gene_names,protein_name",
            "query": query,
        }
    )
    url = f"https://rest.uniprot.org/uniprotkb/stream?{params}"
    req = urllib.request.Request(url, headers={"User-Agent": "OCPT-Ludmila-figure5/1.0"})
    with urlopen_uniprot(req, timeout=120) as resp:
        text = resp.read().decode("utf-8", errors="replace")
    out: dict[str, UniProtMeta] = {}
    reader = csv.reader(io.StringIO(text), delimiter="\t")
    rows = list(reader)
    if len(rows) < 2:
        return out
    header = [h.strip() for h in rows[0]]
    hmap = {h: i for i, h in enumerate(header)}
    i_entry = hmap.get("Entry", 0)
    i_gp = hmap.get("Gene Names (primary)")
    i_gn = hmap.get("Gene Names")
    i_pn = hmap.get("Protein names")
    for row in rows[1:]:
        if i_entry >= len(row):
            continue
        acc = row[i_entry].strip()
        if not acc:
            continue
        gp = row[i_gp].strip() if i_gp is not None and i_gp < len(row) else ""
        gn = row[i_gn].strip() if i_gn is not None and i_gn < len(row) else ""
        pn = row[i_pn].strip() if i_pn is not None and i_pn < len(row) else ""
        out[acc] = (gp, gn, pn)
    return out


def ensure_uniprot_meta(accessions: set[str]) -> dict[str, UniProtMeta]:
    """Fill cache and return acc -> (gene_primary, gene_names, protein_name). Refetches empty rows."""
    cache = load_uniprot_cache()
    missing: list[str] = []
    for a in accessions:
        if not a:
            continue
        meta = cache.get(a)
        if meta is None or not _uniprot_row_nonempty(meta):
            missing.append(a)
    updated = False
    for i in range(0, len(missing), UNIPROT_CHUNK):
        chunk = missing[i : i + UNIPROT_CHUNK]
        try:
            got = fetch_uniprot_chunk(chunk)
        except Exception:
            got = {}
        for acc in chunk:
            cache[acc] = got.get(acc, ("", "", ""))
        updated = True
        if i + UNIPROT_CHUNK < len(missing):
            time.sleep(UNIPROT_SLEEP_S)
    if updated:
        write_uniprot_cache(cache)
    return {a: cache.get(a, ("", "", "")) for a in accessions if a}


def mutated_feat_has_diamond_hit(feat: dict, locus_to_sseqid: dict[str, str]) -> bool:
    """True if this CDS has a DIAMOND best hit (UniProt-style annotation); skip text labels for these."""
    lt = (feat.get("locus_tag") or "").strip()
    g = (feat.get("gene") or "").strip()
    if lt and lt in locus_to_sseqid:
        return True
    if g and g in locus_to_sseqid:
        return True
    return False


def gene_match_key(name: str) -> set[str]:
    keys: set[str] = set()
    n = name.strip()
    if n:
        keys.add(n)
    for part in n.split("/"):
        p = part.strip()
        if p:
            keys.add(p)
    return keys


def mutated_gene_names(rep_rows: dict[int, list[dict]]) -> set[str]:
    out: set[str] = set()
    for rows in rep_rows.values():
        for r in rows:
            g = (r.get("gene_name") or "").strip()
            if g:
                out |= gene_match_key(g)
            lt = (r.get("locus_tag") or "").strip()
            if lt:
                out |= gene_match_key(lt)
    return out


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


def as_int_pos(row: dict) -> int | None:
    try:
        return int(float(row.get("position", "")))
    except Exception:
        return None


def hypermutator_zoom_xlim(
    mut_rows: list[dict], genome_len: float, max_span_frac: float = 0.45
) -> tuple[float, float] | None:
    """Zoom x-axis when point/small variants occupy a minority of the genome (e.g. PedCF).

    Skips megabase-scale large_deletion / large_insertion rows so one whole-genome DEL
    does not force full-length x-axis while SNPs cluster in a smaller window.
    """
    coords: list[int] = []
    for r in mut_rows:
        cat = (r.get("mutation_category") or "").lower()
        if "large_deletion" in cat or "large_insertion" in cat:
            continue
        p = as_int_pos(r)
        if p is not None and 1 <= p <= genome_len:
            coords.append(p)
        ps_raw = (r.get("position_start") or "").strip()
        pe_raw = (r.get("position_end") or "").strip()
        try:
            ps = int(float(ps_raw)) if ps_raw else None
            pe = int(float(pe_raw)) if pe_raw else None
        except (TypeError, ValueError):
            ps, pe = None, None
        if ps is not None and pe is not None and pe >= ps:
            if (pe - ps) <= max_span_frac * genome_len:
                coords.extend([ps, pe])
        elif ps is not None and 1 <= ps <= genome_len:
            coords.append(ps)
    if len(coords) < 2:
        return None
    lo, hi = min(coords), max(coords)
    span = hi - lo
    if span <= 0 or span / genome_len > max_span_frac:
        return None
    pad = max(2500.0, 0.035 * genome_len, 0.1 * span)
    x0 = max(0.0, float(lo) - pad)
    x1 = min(float(genome_len), float(hi) + pad)
    if x1 - x0 < 1e3:
        return None
    return (x0, x1)


def feature_is_mutated(feat: dict, mutated: set[str]) -> bool:
    keys: set[str] = set()
    for field in (feat.get("gene"), feat.get("locus_tag")):
        if not field:
            continue
        keys |= gene_match_key(str(field).strip())
    return bool(keys & mutated)


def intergenic_mutation_positions(
    rep_rows: dict[int, list[dict]], feats: list[dict], genome_len: float
) -> list[int]:
    """Positions with a variant call that lie outside every parsed feature interval."""
    intervals = [(f["start"], f["end"]) for f in feats]

    def inside_any(p: int) -> bool:
        return any(s <= p <= e for s, e in intervals)

    seen: set[int] = set()
    ordered: list[int] = []
    for rows in rep_rows.values():
        for r in rows:
            cat = (r.get("mutation_category") or "").lower()
            if "large_deletion" in cat or "large_insertion" in cat:
                continue
            p = as_int_pos(r)
            if p is None or p < 1 or p > genome_len:
                continue
            if inside_any(p):
                continue
            if p not in seen:
                seen.add(p)
                ordered.append(p)
    return ordered


def main_map_labeled_cds_intervals(
    rep_to_sample: dict[int, str],
    cds_list: list[dict],
    locus_to_sseqid: dict[str, str],
) -> set[tuple[int, int]]:
    """(start, end) of CDS that receive text labels on the main multi-rep map (diamond-skipped excluded)."""
    rep_rows_main: dict[int, list[dict]] = {}
    for rep, sid in rep_to_sample.items():
        if sid in HYPERMUTATOR_EXCLUDE_FROM_MAIN_MAP:
            rep_rows_main[rep] = []
            continue
        rep_rows_main[rep] = load_mutations(sid)
    mut_main = mutated_gene_names(rep_rows_main)
    feats = [f for f in cds_list if feature_is_mutated(f, mut_main)]
    feats_s = sorted(feats, key=lambda f: (f["start"], f["end"]))
    out: set[tuple[int, int]] = set()
    for feat in feats_s:
        if mutated_feat_has_diamond_hit(feat, locus_to_sseqid):
            continue
        out.add((feat["start"], feat["end"]))
    return out


# Gene symbols (after lowercasing and stripping _1 / _2 suffixes) for surface / wall / phage-host factors
_PHAGE_HOST_FACTOR_GENE_BASES = frozenset(
    {
        "tago",
        "tagu",
        "tagh",
        "tagg",
        "tarb",
        "tarf",
        "tark",
        "tarl",
        "taro",
        "atl",
        "spa",
        "clfa",
        "clfb",
        "fnba",
        "fnbb",
        "sdrc",
        "sdrd",
        "sdre",
        "walr",
        "walk",
        "mprf",
        "dlta",
        "dltb",
        "dltc",
        "dltd",
        "ltas",
        "fmta",
        "saer",
        "saes",
        "srap",
        "lyto",
        "lyra",
    }
)


def _gene_base_symbol(g: str | None) -> str:
    if not g:
        return ""
    s = g.strip().lower()
    return re.sub(r"_\d+$", "", s)


def feat_matches_phage_host_factor(feat: dict) -> bool:
    """True if CDS is annotated as a common S. aureus phage receptor / surface / WTA–LTA factor."""
    g = feat.get("gene")
    if g:
        base = _gene_base_symbol(g)
        if base in _PHAGE_HOST_FACTOR_GENE_BASES:
            return True
    prod = (feat.get("product") or "").lower()
    if prod:
        needles = (
            "teichoic",
            "lipoteichoic",
            "autolysin",
            "immunoglobulin g-binding protein a",
            "clumping factor",
            "fibronectin-binding",
            "serine-aspartate repeat-containing",
            "n-acetylglucosaminyl 1-phosphate transferase",
            "d-alanine--d-alanyl carrier protein ligase",
            "teichoic acid d-alanyltransferase",
            "lipoteichoic acid synthase",
            "wall teichoic",
        )
        if any(n in prod for n in needles):
            return True
    return False


def _label_intervals_fit(
    intervals: list[tuple[float, float]], left: float, right: float, gap: float
) -> bool:
    for lo, hi in intervals:
        if not (right + gap <= lo or left - gap >= hi):
            return False
    return True


def layout_mutated_gene_labels(
    x_centers: list[float],
    texts: list[str],
    genome_len: float,
    fig_w_inches: float,
    fontsize_pt: float = GENE_LABEL_FONTSIZE_PT,
    *,
    gap_mult: float = 1.0,
    step_mult: float = 1.0,
) -> list[tuple[float, int]]:
    """Assign each label an x position and row index to reduce overlap (genomic order)."""
    n = len(x_centers)
    if n == 0:
        return []

    def half_width_bp(i: int) -> float:
        ch = max(len(texts[i]), 3)
        axis_in = max(fig_w_inches * 0.82, 1.0)
        char_in = (fontsize_pt / 72.0) * 0.56
        w = (ch * char_in / axis_in) * genome_len * 0.5
        if x_centers[i] > genome_len * RIGHT_PACK_GENOME_FRAC:
            return w * 1.12
        return w

    gap = max(7200.0, genome_len * 0.00175) * gap_mult
    step = max(2200.0, genome_len / 520.0) * step_mult
    order = sorted(range(n), key=lambda i: x_centers[i])
    rows_intervals: list[list[tuple[float, float]]] = []
    x_placed = [0.0] * n
    row_of = [0] * n

    for ii in order:
        xm = x_centers[ii]
        hw = half_width_bp(ii)
        max_c = genome_len - hw
        min_c = hw

        def try_on_row(r: int, x_try: float) -> bool:
            x_c = min(max(x_try, min_c), max_c)
            left, right = x_c - hw, x_c + hw
            if left < -genome_len * 0.001 or right > genome_len * 1.005:
                return False
            while len(rows_intervals) <= r:
                rows_intervals.append([])
            if not _label_intervals_fit(rows_intervals[r], left, right, gap):
                return False
            rows_intervals[r].append((left, right))
            x_placed[ii] = x_c
            row_of[ii] = r
            return True

        placed = False
        max_k = int(genome_len / step) + 16
        if xm > genome_len * RIGHT_PACK_GENOME_FRAC:
            k_trials = [0] + [-(j) for j in range(1, max_k)] + list(range(1, max_k))
        else:
            k_trials = list(range(max_k))
        for r in range(len(rows_intervals) + 1):
            for k in k_trials:
                if try_on_row(r, xm + k * step):
                    placed = True
                    break
            if placed:
                break
        if not placed:
            r = len(rows_intervals)
            while len(rows_intervals) <= r:
                rows_intervals.append([])
            x_c = min(max(xm, min_c), max_c)
            left, right = x_c - hw, x_c + hw
            rows_intervals[r].append((left, right))
            x_placed[ii] = x_c
            row_of[ii] = r

    return [(x_placed[i], row_of[i]) for i in range(n)]


def assign_leader_horizontal_lanes(
    xm_genes: list[float],
    x_labels: list[float],
    row_of: list[int],
    genome_len: float,
) -> list[int]:
    """Give each leader a lane so horizontal segments at the same label row do not overlap in x."""
    n = len(xm_genes)
    if n == 0:
        return []
    x_margin = max(3800.0, genome_len * 0.00012)
    lane_h = [0] * n
    by_row: dict[int, list[int]] = defaultdict(list)
    for i in range(n):
        by_row[row_of[i]].append(i)

    for _r, idxs in by_row.items():
        idxs_sorted = sorted(idxs, key=lambda i: min(xm_genes[i], x_labels[i]))
        occupied: list[list[tuple[float, float]]] = []

        for i in idxs_sorted:
            lo = min(xm_genes[i], x_labels[i])
            hi = max(xm_genes[i], x_labels[i])
            placed = False
            for L in range(len(occupied) + 1):
                while len(occupied) <= L:
                    occupied.append([])
                ok = True
                for olo, ohi in occupied[L]:
                    if not (hi + x_margin <= olo or lo - x_margin >= ohi):
                        ok = False
                        break
                if ok:
                    occupied[L].append((lo, hi))
                    lane_h[i] = L
                    placed = True
                    break
            if not placed:
                L = len(occupied)
                occupied.append([(lo, hi)])
                lane_h[i] = L

    return lane_h


def assign_label_text_stagger(
    layout: list[tuple[float, int]],
    lane_h: list[int],
    n_tiers: int,
    tier_dy: float,
) -> list[float]:
    """Per-label extra downward offset (data coords) so text tiers stagger along the genome."""
    n = len(layout)
    if n == 0:
        return []
    tier_idx = [0] * n
    groups: dict[tuple[int, int], list[int]] = defaultdict(list)
    for i in range(n):
        x_label, row = layout[i]
        groups[(row, lane_h[i])].append(i)
    for _key, idxs in groups.items():
        for j, ii in enumerate(sorted(idxs, key=lambda k: layout[k][0])):
            tier_idx[ii] = j % max(1, n_tiers)
    return [tier_idx[i] * tier_dy for i in range(n)]


def cluster_mutation_x_offsets(
    events: list[tuple[int, str]],
    genome_len: float,
) -> list[tuple[float, str, int]]:
    """Spread nearby genomic positions horizontally (bp) so shaped markers do not overlap."""
    if not events:
        return []
    gap_max = max(3800.0, genome_len * 0.00042)
    step = max(1600.0, genome_len / 8500.0)
    ev = sorted(events, key=lambda e: (e[0], e[1]))
    clusters: list[list[tuple[int, str]]] = []
    cur: list[tuple[int, str]] = [ev[0]]
    for bp, tt in ev[1:]:
        if bp - cur[-1][0] <= gap_max:
            cur.append((bp, tt))
        else:
            clusters.append(cur)
            cur = [(bp, tt)]
    clusters.append(cur)
    out: list[tuple[float, str, int]] = []
    for cl in clusters:
        n = len(cl)
        xs = [p for p, _ in cl]
        mean_x = sum(xs) / n
        cl_step = step
        if mean_x > genome_len * 0.62:
            cl_step = step * 1.55
        elif mean_x > genome_len * RIGHT_PACK_GENOME_FRAC:
            cl_step = step * 1.22
        for i, (bp, tt) in enumerate(cl):
            x = mean_x + (i - (n - 1) / 2.0) * cl_step
            out.append((float(x), tt, bp))
    return out


def save_legend_once() -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    leg_path = OUT_DIR / f"{LEGEND_BASENAME}.png"
    fig, ax = plt.subplots(figsize=(3.2, 2.2))
    ax.set_axis_off()
    items = [
        *[Patch(facecolor=c, edgecolor="#333333", linewidth=0.5, label=t) for t, c in TYPE_COLORS.items()],
        Patch(facecolor=MUTATED_GENE_COLOR, edgecolor="#1e6b3a", linewidth=0.5, label="Mutated gene"),
        Patch(facecolor=TRACK_GREY, edgecolor=TRACK_EDGE, linewidth=0.5, label="Genome backbone"),
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
    fig.savefig(OUT_DIR / f"{LEGEND_BASENAME}.pdf", bbox_inches="tight")
    plt.close(fig)
    return leg_path


def save_shaped_ticks_legend_once() -> Path:
    OUT_DIR_SHAPED.mkdir(parents=True, exist_ok=True)
    leg_path = OUT_DIR_SHAPED / f"{LEGEND_SHAPED_BASENAME}.png"
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    ax.set_axis_off()
    edge = "#2a2a2a"
    items: list[Line2D | Patch] = []
    for t, c in TYPE_COLORS.items():
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
            mk, ms = TYPE_MARKERS.get(t, TYPE_MARKER_DEFAULT)
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
        Patch(facecolor=MUTATED_GENE_COLOR, edgecolor="#1e6b3a", linewidth=0.5, label="Mutated gene")
    )
    items.append(Patch(facecolor=TRACK_GREY, edgecolor=TRACK_EDGE, linewidth=0.5, label="Genome backbone"))
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
    fig.savefig(OUT_DIR_SHAPED / f"{LEGEND_SHAPED_BASENAME}.pdf", bbox_inches="tight")
    plt.close(fig)
    return leg_path


def save_hypermutator_legend_once() -> Path:
    """Standalone legend for hypermutator genome maps (not embedded in strain figures)."""
    OUT_DIR_HYPERMUTATOR.mkdir(parents=True, exist_ok=True)
    leg_path = OUT_DIR_HYPERMUTATOR / f"{LEGEND_HYPERMUTATOR_BASENAME}.png"
    fig, ax = plt.subplots(figsize=(3.2, 2.2))
    ax.set_axis_off()
    items = [
        *[Patch(facecolor=c, edgecolor="#333333", linewidth=0.5, label=t) for t, c in TYPE_COLORS.items()],
        Patch(facecolor=MUTATED_GENE_COLOR, edgecolor="#1e6b3a", linewidth=0.5, label="Mutated gene"),
        Patch(facecolor=TRACK_GREY, edgecolor=TRACK_EDGE, linewidth=0.5, label="Genome backbone"),
    ]
    ax.legend(
        handles=items,
        loc="center",
        frameon=True,
        fancybox=False,
        edgecolor="#cccccc",
        fontsize=mpl.rcParams["axes.labelsize"],
        title="Hypermutator genome maps",
        title_fontsize=mpl.rcParams["axes.titlesize"],
    )
    fig.savefig(leg_path, dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR_HYPERMUTATOR / f"{LEGEND_HYPERMUTATOR_BASENAME}.pdf", bbox_inches="tight")
    plt.close(fig)
    return leg_path


def plot_strain(
    host_id: str,
    rep_to_sample: dict[int, str],
    ref_gbk: Path,
    *,
    shaped_mutation_markers: bool = False,
    hypermutator_only: bool = False,
) -> Path | None:
    if hypermutator_only:
        shaped_mutation_markers = False

    genome_len, cds_list = parse_cds_features(ref_gbk)
    if genome_len <= 0:
        return None

    ref_id = ref_gbk.parent.name
    diamond_tsv = DIAMOND_DIR / ref_id / f"{ref_id}_diamond_best.tsv"
    locus_to_sseqid = load_diamond_qseqid_map(diamond_tsv)

    use_ocpt_clinical = (
        not hypermutator_only
        and host_id.upper().replace("-", "").replace(" ", "") == "OCPT04SA01"
        and OCPT04SA01_CLINICAL_MUT_TSV.is_file()
    )

    rep_hyper: int | None = None
    if hypermutator_only:
        rep_hyper = HYPERMUTATOR_REP_BY_HOST.get(host_id)
        if rep_hyper is None:
            return None
        sid_hm = rep_to_sample.get(rep_hyper)
        if not sid_hm or sid_hm not in HYPERMUTATOR_SAMPLES:
            return None
        rep_rows = {rep_hyper: load_mutations(sid_hm)}
    else:
        rep_rows = {}
        for rep, sid in rep_to_sample.items():
            if sid in HYPERMUTATOR_EXCLUDE_FROM_MAIN_MAP:
                rep_rows[rep] = []
                continue
            rep_rows[rep] = load_mutations(sid)

        if use_ocpt_clinical:
            rep_rows[OCPT04SA01_CLINICAL_REP] = load_mutations(
                OCPT04SA01_CLINICAL_MUT_TSV.stem
            )

    mutated = mutated_gene_names(rep_rows)

    ocpt_hm_post_treatment_intervals: set[tuple[int, int]] | None = None
    if hypermutator_only and host_id == "OCPT04SA01" and rep_hyper is not None:
        ocpt_hm_post_treatment_intervals = main_map_labeled_cds_intervals(
            rep_to_sample, cds_list, locus_to_sseqid
        )

    y_track = 0.0
    track_half = 0.30
    y_reps_full = (
        {1: 6.0, 2: 5.0, 3: 4.0, 4: 3.0, 5: 2.0, 6: 1.0}
        if use_ocpt_clinical
        else {1: 5.0, 2: 4.0, 3: 3.0, 4: 2.0, 5: 1.0}
    )
    y_reps = (
        {rep_hyper: 1.0}
        if hypermutator_only and rep_hyper is not None
        else y_reps_full
    )
    reps_to_plot = (
        [rep_hyper]
        if hypermutator_only and rep_hyper is not None
        else (list(range(1, 7)) if use_ocpt_clinical else list(range(1, 6)))
    )
    tick_half = 0.24 if shaped_mutation_markers else 0.2
    grey_rep4 = host_id.upper().replace("-", "").replace(" ", "") in {
        h.upper().replace("-", "").replace(" ", "") for h in GREY_OUT_REP4_HOSTS
    }
    if hypermutator_only:
        out_dir = OUT_DIR_HYPERMUTATOR
        map_suffix = ""
    else:
        out_dir = OUT_DIR_SHAPED if shaped_mutation_markers else OUT_DIR
        map_suffix = "_shaped_ticks" if shaped_mutation_markers else ""

    mutated_feats = [f for f in cds_list if feature_is_mutated(f, mutated)]

    mutated_sorted = sorted(mutated_feats, key=lambda f: (f["start"], f["end"]))
    label_xm: list[float] = []
    label_texts: list[str] = []
    for feat in mutated_sorted:
        if mutated_feat_has_diamond_hit(feat, locus_to_sseqid):
            continue
        s, e = feat["start"], feat["end"]
        if ocpt_hm_post_treatment_intervals is not None:
            if (s, e) not in ocpt_hm_post_treatment_intervals and not feat_matches_phage_host_factor(
                feat
            ):
                continue
        if (
            host_id.upper().replace("-", "").replace(" ", "") == "OCPT04SA01"
            and not feature_has_named_gene(feat)
        ):
            continue
        label_xm.append((s + e) / 2.0)
        label_texts.append(label_mutated_cds(feat))

    ocpt_hm_dense = hypermutator_only and host_id == "OCPT04SA01"
    ocpt_clinical_w = 1.0
    ocpt_clinical_h = 1.0
    if not hypermutator_only and use_ocpt_clinical:
        n_clinical = len(rep_rows.get(OCPT04SA01_CLINICAL_REP, []))
        if n_clinical > 0:
            # User asked to show every variant even if the canvas is very wide/tall.
            ocpt_clinical_w = 1.10 + min(3.0, 0.072 * float(n_clinical))
            ocpt_clinical_h = 1.04 + min(1.15, 0.034 * float(max(0, n_clinical - 4)))
    label_fs_pt = (
        GENE_LABEL_FONTSIZE_PT * OCPT_HM_LABEL_FONT_SCALE if ocpt_hm_dense else GENE_LABEL_FONTSIZE_PT
    )
    gap_m = OCPT_HM_LAYOUT_GAP_MULT if ocpt_hm_dense else 1.0
    step_m = OCPT_HM_LAYOUT_STEP_MULT if ocpt_hm_dense else 1.0
    fig_w = (
        STRAIN_MAP_FIG_W_IN
        * (OCPT_HM_FIG_W_MULT if ocpt_hm_dense else 1.0)
        * ocpt_clinical_w
    )
    layout = layout_mutated_gene_labels(
        label_xm,
        label_texts,
        float(genome_len),
        fig_w,
        label_fs_pt,
        gap_mult=gap_m,
        step_mult=step_m,
    )
    x_placed = [layout[i][0] for i in range(len(layout))]
    row_placed = [layout[i][1] for i in range(len(layout))]
    lane_h = assign_leader_horizontal_lanes(
        label_xm, x_placed, row_placed, float(genome_len)
    )
    text_stagger = assign_label_text_stagger(
        layout,
        lane_h,
        LABEL_TEXT_STAGGER_TIERS,
        LABEL_TEXT_STAGGER_DY,
    )
    leader_lane_dy = 0.062
    fig_h = STRAIN_MAP_FIG_H_IN * ocpt_clinical_h
    if hypermutator_only:
        fig_h = STRAIN_MAP_FIG_H_HYPERMUTATOR_IN * (
            OCPT_HM_FIG_H_MULT if ocpt_hm_dense else 1.0
        )
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    y_bar_bottom = y_track - track_half
    # Extra stub drop adds padding between genome track / x region and gene labels
    leader_stub_drop = 0.28
    label_row_step = OCPT_HM_LABEL_ROW_STEP if ocpt_hm_dense else 0.52
    text_below_arm = 0.26

    # Single continuous grey genome backbone (no per-gene grey segments)
    ax.barh(
        y_track,
        genome_len,
        left=0,
        height=2 * track_half,
        color=TRACK_GREY,
        edgecolor=TRACK_EDGE,
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
            color=MUTATED_GENE_COLOR,
            edgecolor="#1e6b3a",
            linewidth=0.35,
            zorder=3,
        )

    for p in intergenic_mutation_positions(rep_rows, cds_list, float(genome_len)):
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
            color=MUTATED_GENE_COLOR,
            edgecolor="#1e6b3a",
            linewidth=0.35,
            zorder=3,
        )

    for i in range(len(label_xm)):
        xm_gene = label_xm[i]
        lab = label_texts[i]
        x_label, row = layout[i]
        lh = lane_h[i]
        y_arm = (
            y_bar_bottom
            - leader_stub_drop
            - row * label_row_step
            - lh * leader_lane_dy
        )
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
            fontsize=label_fs_pt,
            color="#222222",
            clip_on=False,
            zorder=4,
        )

    for rep in reps_to_plot:
        y = y_reps[rep]
        sid = rep_to_sample.get(rep, "")
        rep_line = "#cfcfcf" if (grey_rep4 and rep == 4) else "#e8e8e8"
        ax.axhline(y, xmin=0, xmax=1, color=rep_line, linewidth=0.75, zorder=1)

        if not hypermutator_only and grey_rep4 and rep == 4:
            draw_excluded_replicate_tag(ax, float(genome_len), y)

        if sid in HYPERMUTATOR_EXCLUDE_FROM_MAIN_MAP:
            if not hypermutator_only:
                continue

        rows = rep_rows.get(rep, [])
        y_track_top = y_track + track_half
        grey_row = grey_rep4 and rep == 4

        if shaped_mutation_markers:
            events: list[tuple[int, str]] = []
            for r in rows:
                p = as_int_pos(r)
                if p is None or p < 1 or p > genome_len:
                    continue
                t = (r.get("type") or "SNP").strip().upper()
                events.append((p, t))
            for x_draw, t, _bp in cluster_mutation_x_offsets(events, float(genome_len)):
                color = (
                    TYPE_COLORS_MUTED.get(t, "#6b6b6b")
                    if grey_row
                    else TYPE_COLORS.get(t, "#4d4d4d")
                )
                conn_lw, conn_alpha = (1.0, 0.42) if grey_row else (1.05, 0.52)
                if grey_row:
                    tick_lw, tick_alpha = 1.85, 0.78
                else:
                    tick_lw, tick_alpha = 2.0, 0.96
                ax.plot(
                    [x_draw, x_draw],
                    [y_track_top, y - tick_half * 0.9],
                    linestyle=":",
                    color=color,
                    linewidth=conn_lw,
                    alpha=conn_alpha,
                    solid_capstyle="round",
                    zorder=3,
                )
                if t == "SNP":
                    ax.vlines(
                        x_draw,
                        y - tick_half,
                        y + tick_half,
                        color=color,
                        linewidth=tick_lw,
                        alpha=tick_alpha,
                        zorder=4,
                    )
                else:
                    marker, ms = TYPE_MARKERS.get(t, TYPE_MARKER_DEFAULT)
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
            # Single hypermutator row: hundreds of nearby variants stack into one thick line unless spread.
            if hypermutator_only:
                events_hm: list[tuple[int, str]] = []
                for r in rows:
                    p = as_int_pos(r)
                    if p is None or p < 1 or p > genome_len:
                        continue
                    t = (r.get("type") or "SNP").strip().upper()
                    events_hm.append((p, t))
                for x_draw, t, _orig_bp in cluster_mutation_x_offsets(
                    events_hm, float(genome_len)
                ):
                    if grey_row:
                        color = TYPE_COLORS_MUTED.get(t, "#6b6b6b")
                        conn_lw, conn_alpha = 1.0, 0.42
                        tick_lw, tick_alpha = 1.85, 0.78
                    else:
                        color = TYPE_COLORS.get(t, "#4d4d4d")
                        conn_lw, conn_alpha = 1.05, 0.52
                        tick_lw, tick_alpha = 2.0, 0.96
                    y_tick_bottom = y - tick_half
                    ax.plot(
                        [x_draw, x_draw],
                        [y_track_top, y_tick_bottom],
                        linestyle=":",
                        color=color,
                        linewidth=conn_lw,
                        alpha=conn_alpha,
                        solid_capstyle="round",
                        zorder=3,
                    )
                    ax.vlines(
                        x_draw,
                        y - tick_half,
                        y + tick_half,
                        color=color,
                        linewidth=tick_lw,
                        alpha=tick_alpha,
                        zorder=4,
                    )
            else:
                pos_by_type: dict[str, list[int]] = defaultdict(list)
                for r in rows:
                    p = as_int_pos(r)
                    if p is None or p < 1 or p > genome_len:
                        continue
                    t = (r.get("type") or "SNP").strip().upper()
                    pos_by_type[t].append(p)
                for t, positions in pos_by_type.items():
                    if grey_row:
                        color = TYPE_COLORS_MUTED.get(t, "#6b6b6b")
                        conn_lw, conn_alpha = 1.0, 0.42
                        tick_lw, tick_alpha = 1.85, 0.78
                    else:
                        color = TYPE_COLORS.get(t, "#4d4d4d")
                        conn_lw, conn_alpha = 1.05, 0.52
                        tick_lw, tick_alpha = 2.0, 0.96
                    y_tick_bottom = y - tick_half
                    for p in positions:
                        ax.plot(
                            [p, p],
                            [y_track_top, y_tick_bottom],
                            linestyle=":",
                            color=color,
                            linewidth=conn_lw,
                            alpha=conn_alpha,
                            solid_capstyle="round",
                            zorder=3,
                        )
                    ax.vlines(
                        positions,
                        y - tick_half,
                        y + tick_half,
                        color=color,
                        linewidth=tick_lw,
                        alpha=tick_alpha,
                        zorder=4,
                    )

    if label_xm:
        lowest_text_y = y_bar_bottom
        for i in range(len(label_xm)):
            row = layout[i][1]
            lh = lane_h[i]
            y_arm = (
                y_bar_bottom
                - leader_stub_drop
                - row * label_row_step
                - lh * leader_lane_dy
            )
            lowest_text_y = min(
                lowest_text_y, y_arm - text_below_arm - text_stagger[i]
            )
        label_depth = (y_bar_bottom - lowest_text_y) + 0.34
    else:
        label_depth = 0.35
    if shaped_mutation_markers:
        x_pad = max(5000.0, float(genome_len) / 2200.0)
        ax.set_xlim(0, genome_len * 1.002 + 0.4 * x_pad)
    else:
        ax.set_xlim(0, genome_len * 1.002)
    if hypermutator_only and rep_hyper is not None:
        zm = hypermutator_zoom_xlim(
            rep_rows.get(rep_hyper, []), float(genome_len)
        )
        if zm is not None:
            ax.set_xlim(zm[0], zm[1])
        if host_id == "PedCF":
            xl0, xl1 = ax.get_xlim()
            xl1 = min(
                float(genome_len),
                max(xl1, PEDCF_HYPERMUTATOR_XLIM_RIGHT_MIN_MB * 1e6),
            )
            ax.set_xlim(xl0, xl1)
    y_top = (
        1.85
        if hypermutator_only
        else (max(y_reps.values()) + 0.75)
    )
    ax.set_ylim(y_track - track_half - label_depth, y_top)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _p: f"{x / 1e6:.2f}"))
    ax.set_yticks([y_reps[r] for r in reps_to_plot])
    ytick_labels = [
        "Post-Treatment" if (use_ocpt_clinical and r == OCPT04SA01_CLINICAL_REP) else str(r)
        for r in reps_to_plot
    ]
    ax.set_yticklabels(ytick_labels)
    if grey_rep4:
        for i, ytl in enumerate(ax.get_yticklabels()):
            if i == 3:
                ytl.set_color("#6e6e6e")
    ax.tick_params(axis="both", labelsize=mpl.rcParams["xtick.labelsize"])
    ax.tick_params(axis="y", length=0)
    ax.set_title(
        f"S. aureus {host_id}",
        style="italic",
        fontweight="bold",
        pad=14,
    )
    ax.grid(axis="x", color="#f0f0f0", linewidth=0.5)

    if hypermutator_only or host_id in HYPERMUTATOR_REP_BY_HOST:
        draw_hypermutator_axes_badge(ax)

    out_dir.mkdir(parents=True, exist_ok=True)
    out = (
        out_dir / f"{host_id}_hypermutator_genome_mutation_map.png"
        if hypermutator_only
        else out_dir / f"{host_id}_genome_mutation_map{map_suffix}.png"
    )
    fig.tight_layout()
    left_m = 0.14 if use_ocpt_clinical else 0.10
    fig.subplots_adjust(bottom=0.20, left=left_m)
    fig.savefig(out, dpi=300, bbox_inches=None, pad_inches=0)
    fig.savefig(out.with_suffix(".pdf"), bbox_inches=None, pad_inches=0)
    plt.close(fig)
    return out


def main() -> None:
    manuscript_rc()
    for _kind, _d in (("main", OUT_DIR), ("shaped", OUT_DIR_SHAPED),
                      ("hypermutator", OUT_DIR_HYPERMUTATOR)):
        if want_map(_kind):
            _d.mkdir(parents=True, exist_ok=True)
    samples = load_samples()
    groups: dict[str, dict[int, str]] = defaultdict(dict)
    ref_by_host: dict[str, str] = {}

    for sid, row in samples.items():
        if row["sample_type"] == "host_reference":
            ref_by_host[row["host_id"]] = row["sample_id"]
        elif row["sample_type"] == "resistant_mutant":
            rep = (row.get("replicate") or "").strip()
            if rep.isdigit():
                groups[row["host_id"]][int(rep)] = sid

    outputs: list[Path] = []
    if want_map("main"):
        outputs.append(save_legend_once())
    if want_map("shaped"):
        outputs.append(save_shaped_ticks_legend_once())
    if want_map("hypermutator"):
        outputs.append(save_hypermutator_legend_once())

    for host_id, rep_map in sorted(groups.items()):
        if ONLY_STRAINS and host_id not in ONLY_STRAINS:
            continue
        ref_id = ref_by_host.get(host_id)
        if not ref_id:
            continue
        gbk = REF_DIR / ref_id / f"{ref_id}.gbk"
        if not gbk.exists():
            print(f"Reference GenBank not found for {host_id}: {gbk}")
            continue
        if want_map("main"):
            out = plot_strain(host_id, rep_map, gbk)
            if out:
                outputs.append(out)
        if want_map("shaped"):
            out_s = plot_strain(host_id, rep_map, gbk, shaped_mutation_markers=True)
            if out_s:
                outputs.append(out_s)

    for hm_host in sorted(HYPERMUTATOR_REP_BY_HOST):
        if not want_map("hypermutator"):
            break
        if ONLY_STRAINS and hm_host not in ONLY_STRAINS:
            continue
        ref_id = ref_by_host.get(hm_host)
        if not ref_id:
            continue
        gbk = REF_DIR / ref_id / f"{ref_id}.gbk"
        if not gbk.exists():
            continue
        rep_map = groups.get(hm_host, {})
        out_hm = plot_strain(hm_host, rep_map, gbk, hypermutator_only=True)
        if out_hm:
            outputs.append(out_hm)

    print("Generated:")
    for p in outputs:
        print(f"- {p}")


if __name__ == "__main__":
    main()
