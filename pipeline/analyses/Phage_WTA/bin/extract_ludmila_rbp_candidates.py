#!/usr/bin/env python3
"""
Extract candidate RBP / tail proteins from a Prokka phage annotation directory.

Heuristics:
  (1) product text matches tail/RBP keywords (case-insensitive), or
  (2) protein length > threshold AND CDS lies within a window of CDS indices of
      any keyword-matched "anchor" gene (genome order from GFF or TSV).

Example (Ludmila):
    python bin/extract_ludmila_rbp_candidates.py \\
        --prokka-dir data/phage_prokka/Ludmila \\
        --out data/phage_rbp/Ludmila_RBP_candidates.fasta
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

# Product substring hits (case-insensitive). "RBPs" covers some database phrasing.
KEYWORDS = (
    "tail fiber",
    "tail spike",
    "tail protein",
    "baseplate",
    "base plate",
    "receptor-binding",
    "receptor binding",
    "rbps",
)

# Long ORFs near tail loci (aa)
DEFAULT_LENGTH_THRESHOLD_AA = 400
# Include CDS at index ± window from any keyword anchor index
DEFAULT_NEIGHBOR_WINDOW = 3


def find_prokka_files(prokka_dir: Path, prefix: str) -> tuple[Path, Path | None, Path | None]:
    """
    Return (faa_path, tsv_path_or_none, gff_path_or_none).
    Prefer {prefix}.faa / {prefix}.tsv / {prefix}.gff; else single *.faa / *.tsv / *.gff.
    """
    faa = prokka_dir / f"{prefix}.faa"
    if not faa.is_file():
        faas = sorted(prokka_dir.glob("*.faa"))
        if len(faas) == 1:
            faa = faas[0]
        else:
            raise FileNotFoundError(
                f"Expected {prokka_dir / (prefix + '.faa')} or exactly one *.faa in {prokka_dir}"
            )

    tsv = prokka_dir / f"{prefix}.tsv"
    if not tsv.is_file():
        tsvs = sorted(prokka_dir.glob("*.tsv"))
        tsv = tsvs[0] if len(tsvs) == 1 else None

    gff = prokka_dir / f"{prefix}.gff"
    if not gff.is_file():
        gffs = sorted(prokka_dir.glob("*.gff"))
        gff = gffs[0] if len(gffs) == 1 else None

    return faa, tsv, gff


def parse_fasta(faa_path: Path) -> dict[str, tuple[str, str]]:
    """locus_tag -> (header_rest, sequence). Locus = first token after '>'."""
    out: dict[str, tuple[str, str]] = {}
    current: str | None = None
    head_rest: str = ""
    seq: list[str] = []
    for line in faa_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            if current is not None:
                out[current] = (head_rest, "".join(seq))
            rest = line[1:].strip()
            tok = rest.split(None, 1)
            current = tok[0]
            head_rest = tok[1] if len(tok) > 1 else ""
            seq = []
        elif current is not None:
            seq.append(line.strip())
    if current is not None:
        out[current] = (head_rest, "".join(seq))
    return out


def load_tsv_products(tsv_path: Path) -> dict[str, str]:
    """Map locus_tag -> product from Prokka .tsv."""
    by_locus: dict[str, str] = {}
    with tsv_path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            return by_locus
        for row in reader:
            tag = (row.get("locus_tag") or row.get("Locus_tag") or "").strip()
            if not tag:
                continue
            prod = (row.get("product") or row.get("Product") or "").strip()
            by_locus[tag] = prod
    return by_locus


def _gff_attrs(attrs: str) -> dict[str, str]:
    d: dict[str, str] = {}
    for field in attrs.strip().split(";"):
        field = field.strip()
        if "=" in field:
            k, v = field.split("=", 1)
            d[k.strip().lower()] = v.strip()
    return d


def parse_gff_cds_order(gff_path: Path) -> list[str]:
    """Ordered locus_tags for CDS features (by contig, then start)."""
    rows: list[tuple[str, int, str]] = []
    with gff_path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _src, ftype, start_s, _end, _s, _strand, _phase, attrs = parts[:9]
            if ftype != "CDS":
                continue
            try:
                start = int(start_s)
            except ValueError:
                continue
            ad = _gff_attrs(attrs)
            locus = ad.get("locus_tag") or ad.get("name") or ""
            if locus:
                rows.append((seqid, start, locus))
    rows.sort(key=lambda x: (x[0], x[1]))
    return [r[2] for r in rows]


def cds_order_from_tsv(tsv_path: Path) -> list[str]:
    """CDS row order in Prokka .tsv (fallback if no GFF)."""
    order: list[str] = []
    with tsv_path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            ftype = (row.get("ftype") or row.get("Ftype") or "").strip()
            if ftype != "CDS":
                continue
            tag = (row.get("locus_tag") or "").strip()
            if tag:
                order.append(tag)
    return order


def product_matches_keywords(product: str) -> bool:
    pl = product.lower()
    for kw in KEYWORDS:
        if kw in pl:
            return True
    return False


def main() -> int:
    ap = argparse.ArgumentParser(description="Extract candidate RBP/tail proteins from Prokka output.")
    ap.add_argument(
        "--prokka-dir",
        type=Path,
        required=True,
        help="Prokka output directory (e.g. data/phage_prokka/Ludmila).",
    )
    ap.add_argument("--out", type=Path, required=True, help="Output FASTA path.")
    ap.add_argument(
        "--length-threshold-aa",
        type=int,
        default=DEFAULT_LENGTH_THRESHOLD_AA,
        help=f"Min length (aa) for proximity rule (default {DEFAULT_LENGTH_THRESHOLD_AA}).",
    )
    ap.add_argument(
        "--neighbor-window",
        type=int,
        default=DEFAULT_NEIGHBOR_WINDOW,
        help=f"CDS index window around keyword anchors (default {DEFAULT_NEIGHBOR_WINDOW}).",
    )
    args = ap.parse_args()

    prokka_dir = args.prokka_dir.resolve()
    if not prokka_dir.is_dir():
        print(f"ERROR: Not a directory: {prokka_dir}", file=sys.stderr)
        return 1

    phage_id = prokka_dir.name
    prefix = phage_id

    # Prokka usually names proteins {prefix}.faa; adjust if your run used another stem.
    faa_path, tsv_path, gff_path = find_prokka_files(prokka_dir, prefix)
    faa_data = parse_fasta(faa_path)
    if not faa_data:
        print(f"ERROR: No sequences in {faa_path}", file=sys.stderr)
        return 1

    if tsv_path is not None:
        products = load_tsv_products(tsv_path)
        if gff_path and gff_path.is_file():
            order = parse_gff_cds_order(gff_path)
        else:
            order = cds_order_from_tsv(tsv_path)
    else:
        print(
            f"Note: no Prokka .tsv under {prokka_dir}; using product text from FAA headers only.",
            file=sys.stderr,
        )
        products = {
            lt: (hr.strip() or "hypothetical protein") for lt, (hr, _) in faa_data.items()
        }
        order = list(faa_data.keys())

    order = [lt for lt in order if lt in faa_data]
    if not order:
        print("ERROR: Could not determine CDS order (GFF CDS lines, TSV CDS rows, or FAA order).", file=sys.stderr)
        return 1

    index_of: dict[str, int] = {}
    for i, lt in enumerate(order):
        index_of.setdefault(lt, i)

    anchors: set[int] = set()
    for lt in order:
        prod = products.get(lt, "")
        if product_matches_keywords(prod):
            anchors.add(index_of[lt])

    neighbor_idxs: set[int] = set()
    w = args.neighbor_window
    for ai in anchors:
        for j in range(max(0, ai - w), ai + w + 1):
            neighbor_idxs.add(j)

    candidates: set[str] = set()
    for lt in order:
        prod = products.get(lt, "")
        seq_rec = faa_data.get(lt)
        if seq_rec is None:
            continue
        _hr, aa = seq_rec
        aa_len = len(aa.replace("*", ""))
        idx = index_of.get(lt, -1)

        if product_matches_keywords(prod):
            candidates.add(lt)
            continue
        if aa_len > args.length_threshold_aa and idx in neighbor_idxs:
            candidates.add(lt)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as out:
        for lt in order:
            if lt not in candidates:
                continue
            prod = products.get(lt, "hypothetical protein")
            safe_prod = prod.replace("|", "/").replace("\t", " ")
            seq_rec = faa_data[lt]
            _hr, aa = seq_rec
            header = f">{phage_id}|gene={lt}|product={safe_prod}"
            out.write(header + "\n")
            chunk = 60
            for i in range(0, len(aa), chunk):
                out.write(aa[i : i + chunk] + "\n")

    print(f"Wrote {len(candidates)} candidate(s) to {args.out.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
