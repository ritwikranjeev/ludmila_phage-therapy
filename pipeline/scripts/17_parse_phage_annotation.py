#!/usr/bin/env python3
"""Parse Ludmila phage FASTA + phold outputs into data/phage/phage_genes.tsv.

Also flags WTA / podovirus receptor–associated phage proteins and writes
data/phage/phage_receptor_associated.tsv for integration with host WTA panel.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import List

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PHAGE_DIR = PROJECT_ROOT / "data" / "phage"
OUT_TSV = PHAGE_DIR / "phage_genes.tsv"
OUT_RECEPTOR = PHAGE_DIR / "phage_receptor_associated.tsv"
IN_PHOLD = PHAGE_DIR / "phold_per_cds_predictions.tsv"
IN_FALLBACK = PHAGE_DIR / "phold_annotations.tsv"

HEADER = [
    "gene_id",
    "contig",
    "start",
    "end",
    "strand",
    "phrog",
    "phold_function",
    "product",
    "function_category",
    "receptor_assoc_tags",
    "is_receptor_associated_candidate",
]

HEADER_REC = [
    "gene_id",
    "contig",
    "start",
    "end",
    "strand",
    "phrog",
    "phold_function",
    "product",
    "tophit_protein",
    "function_category",
    "receptor_assoc_tags",
]


def receptor_association_tags(product: str, phold_function: str, tophit: str) -> List[str]:
    """Keyword / homology hints for WTA-focused podovirus RBPs and baseplate modules."""
    p = (product or "").lower()
    f = (phold_function or "").lower()
    h = (tophit or "").lower()
    text = f"{p} {f} {h}"
    tags: List[str] = []

    if any(
        k in text
        for k in [
            "phi11",
            "gp45",
            "phi 11",
            "phage 11",
        ]
    ):
        tags.append("phi11_gp45_hint")
    if any(
        k in text
        for k in [
            "twort",
            "sa012",
            "φsa012",
            "orf103",
            "orf105",
            "orf 103",
            "orf 105",
        ]
    ):
        tags.append("twort_sa012_orf_hint")
    if any(
        k in text
        for k in [
            "tail fiber",
            "tail spike",
            "tail-spike",
            "receptor-binding",
            "receptor binding",
            "receptor binding protein",
        ]
    ):
        tags.append("tail_fiber_spike_rbp")
    if any(
        k in text
        for k in [
            "baseplate",
            "base plate",
            "wedge",
            "cell wall-binding",
            "cell wall binding",
            "cell wall binding protein",
            "cwbp",
        ]
    ):
        tags.append("baseplate_or_cell_wall_binding")
    if "glycosyltransferase" in text and any(
        k in text for k in ["tar", "tarP", "wta", "teichoic", "glcnac", "glcnac", "wall"]
    ):
        tags.append("tarP_like_or_wta_glycosyltransferase")
    if "glcnac" in text and any(k in text for k in ["wta", "teichoic", "wall", "receptor", "bind"]):
        tags.append("glcnac_wta_binding_context")
    if "receptor" in text and any(k in text for k in ["tail", "baseplate", "fiber", "spike"]):
        tags.append("explicit_receptor_tail_module")
    # phold often labels RBPs generically as "tail" + "tail protein" / "minor tail protein"
    if (phold_function or "").strip().lower() == "tail":
        tags.append("phold_tail_category")
    if (phold_function or "").strip().lower() == "tail" and any(
        k in p for k in ["minor", "fiber", "spike", "receptor", "binding"]
    ):
        tags.append("tail_minor_or_adhesion_language")

    return sorted(set(tags))


def categorize(product: str, phold_function: str) -> str:
    """Map phold function + product text into coarse modules for host–phage integration."""
    p = (product or "").lower()
    fl = (phold_function or "").lower()
    combo = f"{p} {fl}"

    if any(k in combo for k in ["anti-crispr", "anti crispr", "anti_defense", "acr_"]):
        return "anti_defense"
    if any(
        k in combo
        for k in [
            "tail fiber",
            "tail spike",
            "baseplate",
            "receptor-binding",
            "receptor binding",
            "rbp",
            "receptor",
        ]
    ):
        return "receptor_binding"
    if "lysis" in fl or "holin" in p or "lysin" in p or "endolysin" in p:
        return "lysis"
    if "head and packaging" in fl:
        return "head_packaging"
    if fl.strip() == "tail":
        return "tail_structural"
    if "connector" in fl:
        return "connector"
    if "integration and excision" in fl:
        return "integration_excision"
    if "moron" in fl or "host takeover" in fl:
        return "auxiliary_host_takeover"
    if "transcription regulation" in fl:
        return "transcription_regulation"
    if "dna, rna and nucleotide metabolism" in fl:
        return "dna_rna_metabolism"
    if any(k in combo for k in ["polymerase", "helicase", "primase", "ligase"]):
        return "replication"
    if "unknown function" in fl or not fl.strip():
        return "unknown"
    return "other"


def main() -> None:
    src = IN_PHOLD if IN_PHOLD.exists() else IN_FALLBACK
    if not src.exists():
        raise FileNotFoundError(
            f"Missing phold table. Expected {IN_PHOLD} or {IN_FALLBACK}"
        )

    rows: List[dict] = []
    rec_rows: List[dict] = []

    with src.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for r in reader:
            contig = r.get("contig_id", "").strip()
            gene = r.get("cds_id", "").strip()
            phrog = r.get("phrog", "").strip()
            phold_fn = r.get("function", "").strip()
            product = r.get("product", "").strip()
            tophit = r.get("tophit_protein", "").strip()
            fcat = categorize(product, phold_fn)
            rtags = receptor_association_tags(product, phold_fn, tophit)
            tag_str = "|".join(rtags) if rtags else ""
            is_cand = "1" if rtags else "0"

            row = {
                "gene_id": gene,
                "contig": contig,
                "start": r.get("start", "").strip(),
                "end": r.get("end", "").strip(),
                "strand": r.get("strand", "").strip(),
                "phrog": phrog,
                "phold_function": phold_fn,
                "product": product,
                "function_category": fcat,
                "receptor_assoc_tags": tag_str,
                "is_receptor_associated_candidate": is_cand,
            }
            rows.append(row)
            if rtags:
                rec_rows.append(
                    {
                        "gene_id": gene,
                        "contig": contig,
                        "start": r.get("start", "").strip(),
                        "end": r.get("end", "").strip(),
                        "strand": r.get("strand", "").strip(),
                        "phrog": phrog,
                        "phold_function": phold_fn,
                        "product": product,
                        "tophit_protein": tophit,
                        "function_category": fcat,
                        "receptor_assoc_tags": tag_str,
                    }
                )

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_TSV.open("w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=HEADER, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    with OUT_RECEPTOR.open("w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=HEADER_REC, delimiter="\t")
        w.writeheader()
        w.writerows(rec_rows)

    print(f"Read {src.name}: {len(rows)} CDS rows")
    print(f"Wrote {OUT_TSV}")
    print(f"Wrote {OUT_RECEPTOR} ({len(rec_rows)} receptor-associated candidates)")


if __name__ == "__main__":
    main()
