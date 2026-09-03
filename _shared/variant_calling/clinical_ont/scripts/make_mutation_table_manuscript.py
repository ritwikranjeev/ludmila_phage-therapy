#!/usr/bin/env python3
"""
Build a single manuscript / Supplement mutation table from the frozen variant set.

Primary input: tables/variants_for_paper.tsv
Enrichment: tables/mutations.tsv (SnpEff ANN, HGVS, upstream distance)
Optional supplement: tables/variants_highconf_summary.tsv for focus loci missing from paper table
Gene symbols + products: Prokka GFF (locus_tag -> gene, product)

Adjust inclusion rules in CONFIG below.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

# --- Inclusion / QC (tweak for manuscript scope) ---------------------------------
# Defaults align with relaxed ONT high-conf integration (Medaka primary).
MIN_MEDAKA_QUAL = 20.0
MIN_MEDAKA_DP = 10
MIN_AF = 0.85  # applied to Clair3 AF when present; Medaka AF often unset in TSV
# Focus loci merged from mutations.tsv only may have low Medaka QUAL at complex/indels;
# keep a separate floor so validated sites (e.g. fumC upstream, ebh) still appear.
MIN_QUAL_FOCUS_SUPPLEMENT = 2.5
TIERS_INCLUDE = {"Tier 1", "Tier 2"}
# Tier 3 rows kept only if category (normalized) is rRNA-adjacent or locus is forced below
INCLUDE_TIER3_IF_CATEGORY_SUBSTR = ("rrna", "rrna_adjacent")
# Always attempt to add these loci from highconf/mutations if absent after paper pass
SUPPLEMENT_LOCUS_TAGS = {
    "DEBJFDCL_01664",  # agrA
    "DEBJFDCL_00876",  # yclM
    "DEBJFDCL_01461",  # fumC
    "DEBJFDCL_02530",  # mdrP
    "DEBJFDCL_00981",  # ebh
}
PRODUCT_MAX_LEN = 90
# Narrative category overrides (locus_tag -> display category)
CATEGORY_BY_LOCUS = {
    "DEBJFDCL_01664": "regulator",
    "DEBJFDCL_00876": "metabolism",
    "DEBJFDCL_01461": "metabolism",
    "DEBJFDCL_02530": "envelope/transport",
    "DEBJFDCL_00981": "envelope/transport",
}
# Tier shown in manuscript for loci pulled only from mutations.tsv (no high-conf row)
TIER_BY_LOCUS = {
    "DEBJFDCL_01664": "Tier 1",
    "DEBJFDCL_00876": "Tier 2",
    "DEBJFDCL_01461": "Tier 2",
    "DEBJFDCL_02530": "Tier 2",
    "DEBJFDCL_00981": "Tier 2",
}

MANUSCRIPT_COLS = [
    "Chromosome",
    "Position_bp",
    "Reference",
    "Alternate",
    "Gene",
    "Locus_tag",
    "Product_annotation",
    "Effect",
    "Protein_change",
    "Category",
    "Tier",
    "QUAL_Medaka",
    "DP",
    "AF",
    "Callers",
]

SANITY_LOCUS = {
    "DEBJFDCL_01664": "agrA missense (Arg207His)",
    "DEBJFDCL_00876": "yclM missense (Ser4Cys)",
    "DEBJFDCL_01461": "fumC upstream",
    "DEBJFDCL_02530": "mdrP missense",
    "DEBJFDCL_00981": "ebh",
}


def read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def row_key(r: Dict[str, str]) -> Tuple[str, str, str, str]:
    return (r.get("chrom", ""), r.get("pos", ""), r.get("ref", ""), r.get("alt", ""))


def load_focus_extra(path: Path) -> Set[str]:
    tags = set(SUPPLEMENT_LOCUS_TAGS)
    if not path.exists():
        return tags
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        tags.add(line.split()[0])
    return tags


def load_gff_locus_map(gff_path: Path) -> Dict[str, Dict[str, str]]:
    """locus_tag -> gene symbol, long product string."""
    out: Dict[str, Dict[str, str]] = {}
    if not gff_path.exists():
        return out
    gene_re = re.compile(r"gene=([^;]+)")
    name_re = re.compile(r"Name=([^;]+)")
    locus_re = re.compile(r"locus_tag=([^;]+)")
    prod_re = re.compile(r"product=([^;]+)")
    with gff_path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] not in ("CDS", "gene", "rRNA", "tRNA"):
                continue
            attrs = parts[8]
            mloc = locus_re.search(attrs)
            if not mloc:
                continue
            lt = mloc.group(1).strip()
            sym = ""
            mg = gene_re.search(attrs)
            if mg:
                sym = mg.group(1).strip()
            else:
                mn = name_re.search(attrs)
                if mn:
                    sym = mn.group(1).strip()
            prod = ""
            mp = prod_re.search(attrs)
            if mp:
                prod = mp.group(1).strip().replace("%2C", ",")
            if not sym and prod:
                sym = prod
            out[lt] = {"gene_symbol": sym or lt, "product": prod}
    return out


def parse_ann_record(ann: str) -> Dict[str, str]:
    """Parse first SnpEff ANN field in mutations.annotation (full string may be alt|ANN...)."""
    if not ann or "|" not in ann:
        return {}
    bits = ann.split("|")
    out: Dict[str, str] = {}
    if len(bits) > 1:
        out["effect"] = bits[1] if len(bits) > 1 else ""
        out["impact"] = bits[2] if len(bits) > 2 else ""
    if len(bits) > 9:
        out["hgvs_c"] = bits[9] or ""
    if len(bits) > 10:
        out["hgvs_p"] = bits[10] or ""
    if len(bits) > 15 and (bits[15] or "").strip().isdigit():
        out["distance"] = bits[15].strip()
    return out


def upstream_bp_from_hgvs(hgvs_c: str) -> Optional[int]:
    if not hgvs_c:
        return None
    m = re.search(r"c\.-(\d+)", hgvs_c) or re.search(r"n\.-(\d+)", hgvs_c)
    if m:
        return int(m.group(1))
    return None


def normalize_effect(effect: str) -> str:
    e = (effect or "").lower().strip()
    if "missense" in e:
        return "missense"
    if "synonymous" in e:
        return "synonymous"
    if "upstream" in e:
        return "upstream"
    if "downstream" in e:
        return "downstream"
    if "intergenic" in e:
        return "intergenic"
    if "frameshift" in e:
        return "frameshift"
    if "stop_gained" in e or "nonsense" in e:
        return "nonsense"
    return effect or "other"


def normalize_category(cat: str, locus_tag: str) -> str:
    if locus_tag in CATEGORY_BY_LOCUS:
        return CATEGORY_BY_LOCUS[locus_tag]
    c = (cat or "other").strip().lower().replace(" ", "_")
    if "rrna" in c:
        return "rRNA_adjacent"
    mapping = {
        "regulator": "regulator",
        "envelope": "envelope/transport",
        "envelope_transport": "envelope/transport",
        "surface_access": "envelope/transport",
        "metabolism": "metabolism",
        "other": "other",
    }
    return mapping.get(c, cat or "other")


def tier_allowed(tier: str, category: str, locus_tag: str = "") -> bool:
    t = (tier or "").strip()
    c = normalize_category(category or "", locus_tag or "").lower().replace(" ", "_")
    if t in TIERS_INCLUDE:
        return True
    if t == "Tier 3":
        for sub in INCLUDE_TIER3_IF_CATEGORY_SUBSTR:
            if sub in c:
                return True
    return False


def parse_float(x: str, default: float = 0.0) -> float:
    try:
        return float(x) if x not in ("", ".", None) else default
    except (TypeError, ValueError):
        return default


def effective_af(row: Dict[str, str]) -> float:
    c3 = parse_float(row.get("clair3_af") or row.get("clair3_af".lower()) or "0")
    if c3 > 0:
        return c3
    m = parse_float(row.get("medaka_af") or "0")
    if m > 0:
        return m
    if str(row.get("pass_both_callers", "")).lower() == "true":
        return 1.0
    return 0.0


def passes_qc(
    row: Dict[str, str],
    *,
    from_medaka_table: bool = True,
    supplement_focus: bool = False,
) -> bool:
    """from_medaka_table=False: mutations.tsv rows only have qual, no DP/AF."""
    qual = parse_float(row.get("medaka_qual") or row.get("qual") or "0")
    dp = parse_float(row.get("medaka_dp") or "0")
    qmin = MIN_MEDAKA_QUAL
    if supplement_focus:
        qmin = min(MIN_MEDAKA_QUAL, MIN_QUAL_FOCUS_SUPPLEMENT)
    if qual < qmin:
        return False
    if from_medaka_table and dp > 0 and dp < MIN_MEDAKA_DP:
        return False
    if from_medaka_table:
        af = effective_af(row)
        if af > 0 and af < MIN_AF:
            return False
    return True


def mutation_row_to_variant_shape(m: Dict[str, str]) -> Dict[str, str]:
    """Unify mutations.tsv row to highconf-like keys for build_row."""
    raw_ann = m.get("annotation") or ""
    ann = parse_ann_record(raw_ann)
    lt = (m.get("gene") or "").strip()
    eff_detail = (m.get("effect_detail") or ann.get("effect") or "").lower()
    effect = normalize_effect(eff_detail)
    tier = TIER_BY_LOCUS.get(lt, "Tier 2")
    cat = normalize_category("", lt)
    if "rrna" in raw_ann.lower() or "rrna" in (m.get("effect_detail") or "").lower():
        cat = "rRNA_adjacent"
        tier = "Tier 3"
    return {
        "chrom": m.get("chrom", ""),
        "pos": m.get("pos", ""),
        "ref": m.get("ref", ""),
        "alt": m.get("alt", ""),
        "locus_tag": lt,
        "gene": lt,
        "product": m.get("product") or "",
        "effect": effect,
        "impact": ann.get("impact") or "",
        "category": cat,
        "tier": tier,
        "caller_support": "medaka",
        "pass_both_callers": "False",
        "medaka_qual": str(m.get("qual") or ""),
        "medaka_dp": "",
        "medaka_af": "",
        "clair3_af": "",
    }


def format_callers(row: Dict[str, str]) -> str:
    sup = (row.get("caller_support") or "").strip()
    both = str(row.get("pass_both_callers", "")).lower() == "true"
    if sup == "both" or (both and sup in ("", "both")):
        return "Medaka+Clair3"
    if sup == "medaka_only":
        return "Medaka only (high support)"
    if sup:
        return sup.replace("_", " ").title()
    return "Medaka+Clair3" if both else "Medaka"


def strip_p_prefix(hgvs_p: str) -> str:
    p = (hgvs_p or "").strip()
    if p.startswith("p."):
        return p[2:]
    return p


def protein_change_field(
    effect: str,
    mut_row: Optional[Dict[str, str]],
    gsym: str,
    product_col: str,
) -> str:
    ann = parse_ann_record((mut_row or {}).get("annotation") or "")
    eff = normalize_effect(effect)

    if eff == "upstream":
        dist = None
        if ann.get("hgvs_c"):
            dist = upstream_bp_from_hgvs(ann["hgvs_c"])
        if dist is None and ann.get("distance"):
            try:
                dist = int(ann["distance"])
            except ValueError:
                dist = None
        rna = bool(
            re.search(r"\brrna\b", gsym, re.I)
            or "ribosomal RNA" in gsym
            or "ribosomal_rna" in gsym.lower()
        )
        if dist is not None:
            if rna:
                return f"upstream (-{dist} bp from {gsym})"
            return f"upstream (-{dist} bp from start codon of {gsym})"
        return f"upstream of {gsym}"

    if eff == "synonymous":
        pp = strip_p_prefix(ann.get("hgvs_p", "") or (mut_row or {}).get("product") or product_col)
        return pp if pp else "synonymous"

    # missense / indel in coding: prefer mutations.product (p.X) or ANN
    pp = strip_p_prefix((mut_row or {}).get("product") or "")
    if not pp and ann.get("hgvs_p"):
        pp = strip_p_prefix(ann["hgvs_p"])
    if not pp:
        pc = (product_col or "").strip()
        if pc.startswith("p."):
            pp = strip_p_prefix(pc)
    return pp if pp else "—"


def truncate(s: str, n: int) -> str:
    s = (s or "").strip()
    if len(s) <= n:
        return s
    return s[: n - 1] + "…"


def build_row(
    src: Dict[str, str],
    gff: Dict[str, Dict[str, str]],
    mut_by_key: Dict[Tuple[str, str, str, str], Dict[str, str]],
) -> Dict[str, str]:
    lt = (src.get("locus_tag") or src.get("gene") or "").strip()
    if lt.startswith("DEBJFDCL_") is False and (src.get("gene") or "").startswith("DEBJFDCL_"):
        lt = src["gene"].strip()
    ginfo = gff.get(lt, {})
    gsym = ginfo.get("gene_symbol") or lt
    prod_long = truncate(ginfo.get("product") or "", PRODUCT_MAX_LEN)
    if (
        gsym == prod_long
        and "hypothetical" in (gsym or "").lower()
    ):
        gsym = lt
    key = row_key(src)
    mut_row = mut_by_key.get(key)
    effect = normalize_effect(src.get("effect") or "")
    cat = normalize_category(src.get("category") or "other", lt)

    return {
        "Chromosome": src.get("chrom", ""),
        "Position_bp": str(int(float(src.get("pos") or "0"))),
        "Reference": src.get("ref", ""),
        "Alternate": src.get("alt", ""),
        "Gene": gsym,
        "Locus_tag": lt,
        "Product_annotation": prod_long,
        "Effect": effect,
        "Protein_change": protein_change_field(effect, mut_row, gsym, src.get("product") or ""),
        "Category": cat,
        "Tier": (src.get("tier") or "").strip(),
        "QUAL_Medaka": f"{parse_float(src.get('medaka_qual') or src.get('qual') or '0'):.2f}"
        .rstrip("0")
        .rstrip(".")
        or "0",
        "DP": (
            str(int(parse_float(src.get("medaka_dp") or "0")))
            if parse_float(src.get("medaka_dp") or "0") > 0
            else "—"
        ),
        "AF": (
            f"{effective_af(src):.3f}".rstrip("0").rstrip(".")
            if effective_af(src) > 0
            else "—"
        ),
        "Callers": format_callers(src),
    }


def sort_key(r: Dict[str, str]) -> Tuple:
    tier_rank = {"Tier 1": 0, "Tier 2": 1, "Tier 3": 2}.get(r.get("Tier", ""), 3)
    pos = int(r.get("Position_bp") or 0)
    return (r.get("Chromosome", ""), tier_rank, pos)


def write_tsv(path: Path, rows: List[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MANUSCRIPT_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def write_markdown(path: Path, rows: List[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Mutation table (post vs. pre reference)",
        "",
        "| " + " | ".join(MANUSCRIPT_COLS) + " |",
        "| " + " | ".join("---" for _ in MANUSCRIPT_COLS) + " |",
    ]
    for r in rows:
        cells = []
        for c in MANUSCRIPT_COLS:
            v = (r.get(c) or "").replace("|", "\\|").replace("\n", " ")
            cells.append(v)
        lines.append("| " + " | ".join(cells) + " |")
    path.write_text("\n".join(lines) + "\n")


def tex_escape(s: str) -> str:
    rep = {
        "&": r"\&",
        "%": r"\%",
        "_": r"\_",
        "#": r"\#",
        "$": r"\$",
    }
    out = s.replace("\u2212", "-")  # minus → ASCII for pdfLaTeX compatibility
    for a, b in rep.items():
        out = out.replace(a, b)
    return out


def write_latex(path: Path, rows: List[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    ncols = len(MANUSCRIPT_COLS)
    colspec = "l" + "p{2.1cm}" * (ncols - 1)  # first col short; wrap others lightly
    # simpler: all p columns for wrapping
    colspec = "|" + "|".join(["p{1.35cm}"] + ["p{1.9cm}"] * (ncols - 1)) + "|"
    lines = [
        r"\begin{tabular}{" + colspec + "}",
        r"\hline",
        " & ".join(tex_escape(c).replace("_", r"\_") for c in MANUSCRIPT_COLS) + r" \\",
        r"\hline",
    ]
    for r in rows:
        line = " & ".join(tex_escape(str(r.get(c) or "")) for c in MANUSCRIPT_COLS) + r" \\"
        lines.append(line)
    lines.extend([r"\hline", r"\end{tabular}", ""])
    path.write_text("\n".join(lines))


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--paper", type=Path, default=Path("tables/variants_for_paper.tsv"))
    p.add_argument("--highconf", type=Path, default=Path("tables/variants_highconf_summary.tsv"))
    p.add_argument("--mutations", type=Path, default=Path("tables/mutations.tsv"))
    p.add_argument(
        "--gff",
        type=Path,
        default=Path("annotation/snpeff/data/mrsa_ref_pre/genes.gff"),
    )
    p.add_argument(
        "--focus-tags",
        type=Path,
        default=Path("resources/focus_locus_tags.txt"),
    )
    p.add_argument("--tsv-out", type=Path, default=Path("tables/mutation_table_for_manuscript.tsv"))
    p.add_argument("--md-out", type=Path, default=Path("tables/mutation_table_for_manuscript.md"))
    p.add_argument("--tex-out", type=Path, default=Path("tables/mutation_table_for_manuscript.tex"))
    args = p.parse_args()

    paper = read_tsv(args.paper)
    highconf = read_tsv(args.highconf)
    supplement_tags = load_focus_extra(args.focus_tags)

    print("=== Input: variants_for_paper.tsv ===")
    if paper:
        print("colnames:", list(paper[0].keys()))
        for i, r in enumerate(paper[:10]):
            print(f"row{i}:", dict(r))
    else:
        print("(empty or missing)")

    gff = load_gff_locus_map(args.gff)
    mut_list = read_tsv(args.mutations)
    mut_by_key: Dict[Tuple[str, str, str, str], Dict[str, str]] = {}
    for m in mut_list:
        if str(m.get("present_in_post", "")).lower() not in ("true", "1", "yes"):
            continue
        mut_by_key[row_key(m)] = m

    highconf_by_key = {row_key(r): r for r in highconf}

    selected: Dict[Tuple[str, str, str, str], Dict[str, str]] = {}
    for row in paper:
        if not passes_qc(row, from_medaka_table=True):
            continue
        lt = (row.get("locus_tag") or row.get("gene") or "").strip()
        if not tier_allowed(row.get("tier") or "", row.get("category") or "", lt):
            continue
        selected[row_key(row)] = row

    # Supplement: focus loci from highconf, then from Medaka mutations.tsv if still missing
    for lt in sorted(supplement_tags):
        if any((r.get("locus_tag") or r.get("gene")) == lt for r in selected.values()):
            continue
        for r in highconf:
            if (r.get("locus_tag") or "") != lt:
                continue
            if not passes_qc(r, from_medaka_table=True):
                continue
            if not tier_allowed(r.get("tier") or "", r.get("category") or "", lt):
                if lt not in SUPPLEMENT_LOCUS_TAGS:
                    continue
            selected[row_key(r)] = r
            break
        if any((r.get("locus_tag") or r.get("gene")) == lt for r in selected.values()):
            continue
        for m in mut_list:
            if str(m.get("present_in_post", "")).lower() not in ("true", "1", "yes"):
                continue
            if (m.get("gene") or "").strip() != lt:
                continue
            synth = mutation_row_to_variant_shape(m)
            relaxed_qc = (lt in SUPPLEMENT_LOCUS_TAGS) or (
                synth.get("category") == "rRNA_adjacent"
            )
            if not passes_qc(
                synth,
                from_medaka_table=False,
                supplement_focus=relaxed_qc,
            ):
                continue
            selected[row_key(synth)] = synth
            break

    # rRNA-adjacent post variants (e.g. upstream of rRNA) not already included
    for m in mut_list:
        if str(m.get("present_in_post", "")).lower() not in ("true", "1", "yes"):
            continue
        if "rrna" not in (m.get("annotation") or "").lower():
            continue
        synth = mutation_row_to_variant_shape(m)
        if synth.get("category") != "rRNA_adjacent":
            continue
        k = row_key(synth)
        if k in selected:
            continue
        if not passes_qc(synth, from_medaka_table=False, supplement_focus=True):
            continue
        if not tier_allowed(synth["tier"], synth["category"], synth["locus_tag"]):
            continue
        selected[k] = synth

    final_src = list(selected.values())
    out_rows = [build_row(r, gff, mut_by_key) for r in final_src]
    out_rows.sort(key=sort_key)

    print()
    print("=== Manuscript table (first 20 lines) ===")
    for i, r in enumerate(out_rows[:20]):
        print(i, r)

    print()
    print("=== Sanity: expected loci ===")
    have_lt = {r["Locus_tag"] for r in out_rows}
    for lt, desc in SANITY_LOCUS.items():
        if lt in have_lt:
            hit = next(x for x in out_rows if x["Locus_tag"] == lt)
            print(f"OK  {lt} ({desc}): {hit.get('Protein_change', '')}")
        else:
            print(f"MISSING {lt} ({desc}) — check highconf / paper / QC thresholds")

    write_tsv(args.tsv_out, out_rows)
    write_markdown(args.md_out, out_rows)
    write_latex(args.tex_out, out_rows)
    print()
    print("Wrote:", args.tsv_out, args.md_out, args.tex_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
