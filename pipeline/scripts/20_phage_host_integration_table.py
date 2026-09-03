#!/usr/bin/env python3
"""Build a small crosswalk: Ludmila tail/receptor-tagged CDS vs host WTA panel + GWAS hints."""

import csv
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PHAGE_REC = PROJECT_ROOT / "data" / "phage" / "phage_receptor_associated.tsv"
PANEL = PROJECT_ROOT / "config" / "host_receptor_gene_panel.tsv"
HOST_GWAS = (
    PROJECT_ROOT
    / "outputs"
    / "gwas"
    / "post_gwas_fdr"
    / "hostrec_binary__pyseer_host_receptors_binary_fdr.tsv"
)
OUT_DIR = PROJECT_ROOT / "outputs" / "integration"
OUT_CROSS = OUT_DIR / "phage_host_receptor_crosswalk.tsv"
OUT_PANEL = OUT_DIR / "host_receptor_panel_reference.tsv"

# Free-text cue for manuscript notes (not a statistical join)
TAG_TO_HOST_NOTE = {
    "phold_tail_category": "Compare to host WTA loci (tagO/tar*) and pyseer host_rec hits on those genes.",
    "tail_minor_or_adhesion_language": "Minor tail / adhesion language — candidate RBP diversity vs tarS/M/P.",
    "tail_fiber_spike_rbp": "Classic RBP/fiber module — strongest mechanistic prior vs WTA.",
    "phi11_gp45_hint": "ϕ11/gp45-class hint — literature link to WTA/glycan receptors.",
    "twort_sa012_orf_hint": "Twort/ΦSA012-style tail gene language.",
    "baseplate_or_cell_wall_binding": "Baseplate / CW-binding — check wall teichoic acid pathways.",
    "tarP_like_or_wta_glycosyltransferase": "Phage-encoded GT language — compare host tarP/tarS.",
    "glcnac_wta_binding_context": "GlcNAc/WTA binding context.",
    "explicit_receptor_tail_module": "Explicit receptor + tail module wording.",
}


def load_panel_by_tier():
    by_tier = {}
    with PANEL.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            tier = row.get("tier", "").strip()
            cid = row.get("canonical_id", "").strip()
            if not cid:
                continue
            by_tier.setdefault(tier, []).append(cid)
    return by_tier


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    by_tier = load_panel_by_tier()
    tier_summary = "; ".join(f"{k}:[{','.join(v)}]" for k, v in sorted(by_tier.items()))

    with OUT_PANEL.open("w", newline="") as handle:
        w = csv.DictWriter(
            handle,
            fieldnames=["canonical_id", "tier", "synonyms"],
            delimiter="\t",
        )
        w.writeheader()
        with PANEL.open(newline="") as h2:
            for row in csv.DictReader(h2, delimiter="\t"):
                w.writerow(row)
    print(f"Wrote {OUT_PANEL}")

    cross_fields = [
        "phage_gene_id",
        "phold_function",
        "product",
        "function_category",
        "receptor_assoc_tags",
        "integration_note",
        "host_panel_tiers_summary",
    ]
    rows_out = []
    if PHAGE_REC.is_file():
        with PHAGE_REC.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                tags = (row.get("receptor_assoc_tags") or "").split("|")
                notes = []
                for t in tags:
                    t = t.strip()
                    if t in TAG_TO_HOST_NOTE:
                        notes.append(TAG_TO_HOST_NOTE[t])
                rows_out.append(
                    {
                        "phage_gene_id": row.get("gene_id", ""),
                        "phold_function": row.get("phold_function", ""),
                        "product": row.get("product", ""),
                        "function_category": row.get("function_category", ""),
                        "receptor_assoc_tags": row.get("receptor_assoc_tags", ""),
                        "integration_note": " | ".join(notes) if notes else "Review vs host WTA panel + host GWAS.",
                        "host_panel_tiers_summary": tier_summary,
                    }
                )
    with OUT_CROSS.open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=cross_fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows_out)
    print(f"Wrote {OUT_CROSS} ({len(rows_out)} phage CDS)")

    if HOST_GWAS.is_file():
        out_g = OUT_DIR / "host_receptor_gwas_top_binary.tsv"
        with HOST_GWAS.open(newline="") as inp, out_g.open("w", newline="") as outp:
            r = csv.DictReader(inp, delimiter="\t")
            rows = list(r)[:30]
            if r.fieldnames:
                w = csv.DictWriter(outp, fieldnames=r.fieldnames, delimiter="\t")
                w.writeheader()
                w.writerows(rows)
        print(f"Wrote {out_g} (first 30 rows)")
    else:
        print(f"Optional: add {HOST_GWAS.name} for host GWAS snippet (rsync post_gwas_fdr from HPC).")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
