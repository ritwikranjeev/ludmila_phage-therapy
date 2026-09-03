#!/usr/bin/env python3
"""Functional annotation of top Stk2-linked Panaroo genes (|φ| from linkage table).

Requires representative protein sequences from Panaroo outputs:
  - outputs/pangenome/gene_data.csv  (preferred; has protein_sequence + annotation)
  - or outputs/pangenome/combined_protein_CDS.fasta  (headers = Panaroo gene IDs)

Then runs local DIAMOND or BLASTP if installed and a protein DB path is given.

Outputs:
  results/stk2_linked_top_genes.faa
  results/reports/stk2_linked_panaroo_annotations.tsv
  results/reports/stk2_linked_panaroo_annotation_summary.txt
  results/reports/panaroo_display_labels_from_blast.tsv  (merge into config/panaroo_display_overrides.tsv)

Environment / CLI:
  --gene-data-csv PATH
  --combined-protein-fasta PATH
  --diamond-db PATH.dmnd   (runs: diamond blastp)
  --blast-db PATH          (runs: blastp -db)
  --top-n 15

Conservative hit reporting: weak hits flagged; categories heuristic from stitle keywords.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))


def project_root() -> Path:
    return _SCRIPT_DIR.parent


def read_top_genes(link_path: Path, top_n: int, min_c_neg: int = 0) -> pd.DataFrame:
    df = pd.read_csv(link_path, sep="\t")
    df["abs_phi"] = df["phi"].astype(float).abs()
    if "c_stk2neg_genepos" in df.columns and min_c_neg > 0:
        df = df[df["c_stk2neg_genepos"].astype(int) >= min_c_neg]
    df = df.sort_values("abs_phi", ascending=False).drop_duplicates(subset=["gene"]).head(top_n)
    return df


def read_fasta(path: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    name: Optional[str] = None
    chunks: List[str] = []
    with path.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(chunks).replace("*", "")
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            out[name] = "".join(chunks).replace("*", "")
    return out


def load_aa_from_gene_data(csv_path: Path, gene_ids: List[str]) -> Dict[str, Tuple[str, str, str]]:
    """Return gene_id -> (aa_sequence, panaroo_annotation, row_gene_name)."""
    gdf = pd.read_csv(csv_path, low_memory=False)
    # Panaroo column name variants
    name_cols = [c for c in gdf.columns if c in ("gene_name", "Gene", "gff_feature")]
    seq_cols = [
        c
        for c in gdf.columns
        if c.lower()
        in ("protein_sequence", "aa_sequence", "protein", "translation", "prot_sequence")
    ]
    ann_cols = [c for c in gdf.columns if c.lower() in ("annotation", "product", "description")]
    if not ann_cols:
        ann_cols = [c for c in gdf.columns if c.lower() == "description"]
    if not name_cols or not seq_cols:
        raise ValueError(
            f"{csv_path}: need gene name column ({name_cols}) and protein sequence column ({seq_cols}). "
            f"Columns seen: {list(gdf.columns)[:30]}..."
        )
    nk, sk = name_cols[0], seq_cols[0]
    ak = ann_cols[0] if ann_cols else None
    want = set(gene_ids)
    out: Dict[str, Tuple[str, str, str]] = {}
    for r in gdf.itertuples(index=False):
        gid = str(getattr(r, nk)).strip()
        if gid not in want:
            continue
        seq = str(getattr(r, sk)).strip()
        if not seq or seq.lower() in ("nan", "none"):
            continue
        ann = str(getattr(r, ak)).strip() if ak else ""
        out[gid] = (seq.replace("*", ""), ann, gid)
    return out


def load_aa_from_combined_fasta(fasta_path: Path, gene_ids: List[str]) -> Dict[str, Tuple[str, str, str]]:
    fa = read_fasta(fasta_path)
    out = {}
    for g in gene_ids:
        if g in fa:
            out[g] = (fa[g], "", g)
    return out


def write_query_fasta(
    root: Path, gene_rows: pd.DataFrame, seq_map: Dict[str, Tuple[str, str, str]]
) -> Path:
    out_faa = root / "results" / "stk2_linked_top_genes.faa"
    out_faa.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for r in gene_rows.itertuples(index=False):
        gid = str(r.gene)
        if gid not in seq_map:
            continue
        aa, ann, _ = seq_map[gid]
        phi = float(r.phi)
        desc = ann[:120].replace("\n", " ") if ann else "no_panaroo_annotation"
        lines.append(f">{gid} |phi|={abs(phi):.3f} {desc}\n")
        for i in range(0, len(aa), 60):
            lines.append(aa[i : i + 60] + "\n")
    if not lines:
        raise RuntimeError("No sequences written; check gene IDs vs gene_data / combined FASTA.")
    out_faa.write_text("".join(lines), encoding="utf-8")
    return out_faa


def run_diamond(query: Path, db: Path, out_tsv: Path, threads: int = 4) -> None:
    cmd = [
        "diamond",
        "blastp",
        "-d",
        str(db),
        "-q",
        str(query),
        "-o",
        str(out_tsv),
        "--outfmt",
        "6",
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "stitle",
        "--threads",
        str(threads),
    ]
    subprocess.run(cmd, check=True)
    print(f"DIAMOND finished -> {out_tsv}")


def run_blastp(query: Path, db: Path, out_tsv: Path, threads: int = 4) -> None:
    cmd = [
        "blastp",
        "-query",
        str(query),
        "-db",
        str(db),
        "-out",
        str(out_tsv),
        "-outfmt",
        "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle",
        "-num_threads",
        str(threads),
    ]
    subprocess.run(cmd, check=True)
    print(f"blastp finished -> {out_tsv}")


def _remote_blast_stderr_indicates_failure(stderr: str) -> bool:
    if not (stderr or "").strip():
        return False
    low = stderr.lower()
    if "sigxcpu" in low or "internal_error" in low:
        return True
    return "severe error" in low and "search failed" in low


def run_blastp_remote_nr(
    query: Path,
    out_tsv: Path,
    max_target_seqs: int = 1,
    sleep_between: float = 6.0,
    retries: int = 4,
) -> None:
    """NCBI cloud BLAST (no local nr). One query sequence per request — multi-query
    remote jobs often exceed NCBI server CPU limits (SIGXCPU) and return empty output.
    """
    fa = read_fasta(query)
    if not fa:
        raise RuntimeError(f"No sequences in {query}")

    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"
    collected: List[str] = []
    tmp_dir = out_tsv.parent
    tmp_dir.mkdir(parents=True, exist_ok=True)

    for i, (qid, seq) in enumerate(fa.items()):
        if i > 0 and sleep_between > 0:
            time.sleep(sleep_between)
        qpath = tmp_dir / f".remote_blast_q_{i}.faa"
        opath = tmp_dir / f".remote_blast_o_{i}.tsv"
        qpath.write_text(
            f">{qid}\n" + "\n".join(seq[j : j + 60] for j in range(0, len(seq), 60)) + "\n",
            encoding="utf-8",
        )
        if opath.is_file():
            opath.unlink()
        for attempt in range(retries):
            cmd = [
                "blastp",
                "-query",
                str(qpath),
                "-remote",
                "-db",
                "nr",
                "-out",
                str(opath),
                "-outfmt",
                outfmt,
                "-max_target_seqs",
                str(max_target_seqs),
            ]
            proc = subprocess.run(cmd, capture_output=True, text=True)
            if proc.stderr:
                print(proc.stderr, file=sys.stderr, end="")
            size_ok = opath.is_file() and opath.stat().st_size > 0
            bad = proc.returncode != 0 or _remote_blast_stderr_indicates_failure(proc.stderr or "")
            if size_ok and not bad:
                collected.extend(
                    opath.read_text(encoding="utf-8", errors="replace").strip().splitlines()
                )
                break
            if attempt + 1 < retries:
                time.sleep(12.0 * (attempt + 1))
        else:
            print(
                f"[warn] remote blastp: no usable tabular output for {qid} after {retries} tries",
                file=sys.stderr,
            )
        qpath.unlink(missing_ok=True)
        opath.unlink(missing_ok=True)

    out_tsv.write_text("\n".join(collected) + ("\n" if collected else ""), encoding="utf-8")
    if not collected:
        raise RuntimeError(
            "blastp -remote produced no tabular hits for any query; use local DIAMOND/nr or retry later."
        )
    print(f"blastp -remote nr finished -> {out_tsv} ({len(collected)} hit lines)")


def parse_blast_tab(path: Path) -> pd.DataFrame:
    cols = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "stitle",
    ]
    df = pd.read_csv(path, sep="\t", names=cols, comment="#")
    df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")
    df["pident"] = pd.to_numeric(df["pident"], errors="coerce")
    df["bitscore"] = pd.to_numeric(df["bitscore"], errors="coerce")
    return df


def organism_from_stitle(stitle: str) -> str:
    """Extract text in last [...] from NCBI-style subject titles, if present."""
    if "[" not in stitle or "]" not in stitle:
        return ""
    return stitle.rsplit("[", 1)[-1].split("]")[0].strip()[:120]


def classify_hit(stitle: str, pident: float, evalue: float, alen: int) -> Tuple[str, str, str]:
    """Returns (category, hit_confidence, short_reason)."""
    t = (stitle or "").lower()
    weak = []
    if evalue > 1e-3 or pd.isna(evalue):
        weak.append("high_evalue")
    if pident < 35 or pd.isna(pident):
        weak.append("low_pident")
    if alen < 40:
        weak.append("short_aln")
    hit_conf = "low" if len(weak) >= 2 else ("medium" if weak else "high")

    defense_kw = (
        "restriction",
        "crispr",
        "cas ",
        "cas9",
        "toxin",
        "antitoxin",
        "abi",
        "phage",
        "defense",
        "mazg",
        "hok",
        "retron",
        "pd-t",
        "pycsar",
        "gabija",
        "stk2",
        "stk ",
    )
    mobile_kw = (
        "transposase",
        "integrase",
        "recombinase",
        "conjug",
        "mobil",
        "plasmid",
        "is110",
        "is256",
        "resolvase",
    )
    house_kw = (
        "ribosomal",
        "rps",
        "rpl",
        "dna polymerase",
        "rna polymerase",
        "atp synthase",
        "chaperone",
        "elongation factor",
        "trna",
        "gyrase",
    )

    if any(k in t for k in defense_kw):
        cat = "defense_related"
    elif any(k in t for k in mobile_kw):
        cat = "mobile_element"
    elif any(k in t for k in house_kw):
        cat = "housekeeping"
    elif "hypothetical" in t or "uncharacterized" in t or "unknown" in t or "duf" in t:
        cat = "hypothetical"
    else:
        cat = "other"

    if hit_conf == "low":
        cat = "weak_hit"
    return cat, hit_conf, ";".join(weak) if weak else "ok"


def best_hit_per_query(blast_df: pd.DataFrame) -> pd.DataFrame:
    """Keep best hit per qseqid by bitscore then evalue."""
    blast_df = blast_df.copy()
    blast_df["qbase"] = blast_df["qseqid"].astype(str).str.split().str[0]
    blast_df = blast_df.sort_values(["qbase", "bitscore", "evalue"], ascending=[True, False, True])
    return blast_df.groupby("qbase", as_index=False).first()


def build_annotation_table(
    gene_rows: pd.DataFrame,
    seq_map: Dict[str, Tuple[str, str, str]],
    blast_best: Optional[pd.DataFrame],
    tool: str,
) -> pd.DataFrame:
    rows = []
    for r in gene_rows.itertuples(index=False):
        gid = str(r.gene)
        phi = float(r.phi)
        if gid in seq_map:
            aa, pan_ann, _ = seq_map[gid]
            src = "gene_data.csv" if pan_ann else "combined_protein_CDS.fasta"
        else:
            aa, pan_ann, src = "", "", "missing"

        inf_cat = "no_sequence" if not aa else "not_run"
        if aa and tool == "none" and pan_ann:
            pc, _, _ = classify_hit(pan_ann, 100.0, 0.0, max(len(aa), 50))
            if pc != "weak_hit":
                inf_cat = f"panaroo_annotation_only:{pc}"

        row = {
            "panaroo_id": gid,
            "phi": phi,
            "abs_phi": abs(phi),
            "sequence_source": src,
            "aa_length": len(aa) if aa else 0,
            "panaroo_csv_annotation": pan_ann,
            "search_tool": tool,
            "sseqid": "",
            "stitle": "",
            "organism": "",
            "pident": float("nan"),
            "aln_length": float("nan"),
            "evalue": float("nan"),
            "bitscore": float("nan"),
            "inferred_category": inf_cat,
            "hit_confidence": "none",
            "classification_notes": "",
            "display_label_short": "",
        }

        if blast_best is not None and not blast_best.empty:
            hit = blast_best.loc[blast_best["qbase"] == gid]
            if not hit.empty:
                h = hit.iloc[0]
                st = str(h["stitle"])
                pid = float(h["pident"])
                ev = float(h["evalue"])
                alen = int(h["length"])
                cat, hconf, note = classify_hit(st, pid, ev, alen)
                row.update(
                    {
                        "sseqid": str(h["sseqid"]),
                        "stitle": st[:500],
                        "organism": organism_from_stitle(st),
                        "pident": pid,
                        "aln_length": alen,
                        "evalue": ev,
                        "bitscore": float(h["bitscore"]),
                        "inferred_category": cat,
                        "hit_confidence": hconf,
                        "classification_notes": note,
                    }
                )
                # Short figure label: conservative
                if hconf == "high" and cat not in ("weak_hit", "hypothetical"):
                    short = st.split(" >")[0][:70]
                    if len(short) > 68:
                        short = short[:65] + "…"
                    row["display_label_short"] = f"{gid}: {short}"
                elif hconf == "medium":
                    row["display_label_short"] = f"{gid}: putative ({cat.replace('_', ' ')})"
                else:
                    row["display_label_short"] = f"{gid} (weak similarity)"

        if not row["display_label_short"]:
            if pan_ann and len(pan_ann) > 3:
                pa = pan_ann[:60] + ("…" if len(pan_ann) > 60 else "")
                row["display_label_short"] = f"{gid} ({pa})"
            else:
                row["display_label_short"] = f"{gid} (unannotated)"

        rows.append(row)
    return pd.DataFrame(rows)


def write_label_merge_file(ann: pd.DataFrame, path: Path) -> None:
    """Rows for panaroo_display_overrides.tsv (only medium+ blast confidence)."""
    sub = ann[
        (ann["hit_confidence"].isin(["high", "medium"]))
        & (ann["inferred_category"] != "weak_hit")
    ].copy()
    if sub.empty:
        path.write_text(
            "panaroo_id\tdisplay_label\tsource\tconfidence\tnotes\n",
            encoding="utf-8",
        )
        return
    out = sub[
        [
            "panaroo_id",
            "display_label_short",
        ]
    ].rename(
        columns={
            "display_label_short": "display_label",
        }
    )
    out.insert(2, "source", "annotate_stk2_linked_panaroo.py+local_search")
    out["confidence"] = sub["hit_confidence"].values
    out["notes"] = sub["stitle"].str[:200].values
    out.to_csv(path, sep="\t", index=False)


def write_summary_text(ann: pd.DataFrame, path: Path, had_blast: bool) -> None:
    lines = [
        "Stk2-linked Panaroo genes — annotation summary\n",
        "===============================================\n",
        f"Genes in table: {len(ann)}\n",
        f"Local DIAMOND/BLAST run: {'yes' if had_blast else 'no'}\n\n",
        "Category counts (inferred_category):\n",
    ]
    for k, v in ann["inferred_category"].value_counts().items():
        lines.append(f"  {k}: {v}\n")
    lines.append("\nInterpretation (conservative):\n")

    cats = set(ann["inferred_category"])
    defenseish = bool(
        cats.intersection({"defense_related"})
        or ann["stitle"].str.lower().str.contains("recombinase|transposase|integrase", regex=True).any()
    )
    mobileish = bool("mobile_element" in cats or ann["stitle"].str.lower().str.contains("transposase|integrase|plasmid", regex=True).any())

    if not had_blast:
        lines.append(
            "- Search was not run or produced no hits; interpret Panaroo CSV annotations only.\n"
            "- For stronger claims, run DIAMOND/BLASTP with a bacterial protein database.\n"
        )
    if defenseish or mobileish:
        lines.append(
            "- Several linked features match defense- or mobile-element-associated proteins "
            "in reference databases (keyword/heuristic classification; verify top hits manually).\n"
        )
    else:
        lines.append(
            "- No strong automated signal of a single ‘defense island’ type from BLAST titles alone; "
            "many hits may remain hypothetical or weak.\n"
        )
    lines.append(
        "\nStk2 (DefenseFinder) is annotated separately; this table refers to co-segregating Panaroo features.\n"
        "A coherent defense locus is plausible when Stk2+ strains share kdp/xerC-like neighbors "
        "(see linkage φ and optional locus cartoon), but causality is not established here.\n"
    )
    path.write_text("".join(lines), encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project-root", type=Path, default=None)
    ap.add_argument("--linkage-tsv", type=Path, default=None)
    ap.add_argument("--top-n", type=int, default=15)
    ap.add_argument("--min-c-neg", type=int, default=0, help="Min Stk2−&gene+ count (0 = ignore)")
    ap.add_argument("--gene-data-csv", type=Path, default=None)
    ap.add_argument("--combined-protein-fasta", type=Path, default=None)
    ap.add_argument("--diamond-db", type=Path, default=None)
    ap.add_argument("--blast-db", type=Path, default=None)
    ap.add_argument(
        "--blast-remote-nr",
        action="store_true",
        help="Run blastp -remote -db nr (requires network; no local DB).",
    )
    ap.add_argument(
        "--remote-blast-sleep",
        type=float,
        default=6.0,
        help="Seconds between per-sequence NCBI cloud blastp calls (default 6).",
    )
    ap.add_argument(
        "--remote-blast-retries",
        type=int,
        default=4,
        help="Retries per sequence when remote blastp errors or returns no tabular output.",
    )
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument(
        "--raw-hits-tsv",
        type=Path,
        default=None,
        help="Use existing DIAMOND/BLAST tab output (skip search). Same outfmt 6 as produced by this script.",
    )
    args = ap.parse_args()
    root = (args.project_root or project_root()).resolve()
    link = args.linkage_tsv or (root / "results" / "stk2_panaroo_gene_linkage.tsv")
    if not link.is_file():
        raise FileNotFoundError(link)

    top = read_top_genes(link, args.top_n, min_c_neg=args.min_c_neg)
    genes = top["gene"].astype(str).tolist()
    print(f"Top {len(genes)} genes by |φ| (min_c_neg={args.min_c_neg}).", flush=True)

    gd_path = args.gene_data_csv or (root / "outputs" / "pangenome" / "gene_data.csv")
    fa_path = args.combined_protein_fasta or (root / "outputs" / "pangenome" / "combined_protein_CDS.fasta")

    seq_map: Dict[str, Tuple[str, str, str]] = {}
    if gd_path.is_file():
        seq_map = load_aa_from_gene_data(gd_path, genes)
        print(f"Loaded {len(seq_map)} / {len(genes)} AA sequences from {gd_path}")
    if len(seq_map) < len(genes) and fa_path.is_file():
        extra = load_aa_from_combined_fasta(fa_path, genes)
        for k, v in extra.items():
            if k not in seq_map:
                seq_map[k] = v
        print(f"After {fa_path.name}: {len(seq_map)} / {len(genes)} sequences.")

    rep_dir = root / "results" / "reports"
    rep_dir.mkdir(parents=True, exist_ok=True)

    if len(seq_map) == 0:
        msg = (
            f"No sequences found. Place Panaroo `gene_data.csv` or `combined_protein_CDS.fasta` under "
            f"{root / 'outputs' / 'pangenome'} or pass --gene-data-csv / --combined-protein-fasta.\n"
            f"Wrote gene list only -> {rep_dir / 'stk2_linked_top_genes_requested.tsv'}"
        )
        top.to_csv(rep_dir / "stk2_linked_top_genes_requested.tsv", sep="\t", index=False)
        print(msg, file=sys.stderr)
        sys.exit(2)

    query_faa = write_query_fasta(root, top, seq_map)
    print(f"Wrote {query_faa}")

    blast_raw = rep_dir / "stk2_linked_diamond_raw.tsv"
    tool = "none"
    blast_best: Optional[pd.DataFrame] = None
    if args.raw_hits_tsv and args.raw_hits_tsv.is_file():
        blast_best = best_hit_per_query(parse_blast_tab(args.raw_hits_tsv))
        tool = "imported_hits"
        print(f"Loaded existing hits from {args.raw_hits_tsv}", flush=True)
    elif args.diamond_db and args.diamond_db.is_file():
        if not shutil.which("diamond"):
            print("diamond not in PATH; skipping search.", file=sys.stderr)
        else:
            run_diamond(query_faa, args.diamond_db, blast_raw, threads=args.threads)
            blast_best = best_hit_per_query(parse_blast_tab(blast_raw))
            tool = "diamond_blastp"
    elif args.blast_db:
        if not shutil.which("blastp"):
            print("blastp not in PATH; skipping search.", file=sys.stderr)
        else:
            run_blastp(query_faa, args.blast_db, blast_raw, threads=args.threads)
            blast_best = best_hit_per_query(parse_blast_tab(blast_raw))
            tool = "blastp"
    elif args.blast_remote_nr:
        if not shutil.which("blastp"):
            print("blastp not in PATH; skipping remote search.", file=sys.stderr)
        else:
            run_blastp_remote_nr(
                query_faa,
                blast_raw,
                sleep_between=args.remote_blast_sleep,
                retries=args.remote_blast_retries,
            )
            blast_best = best_hit_per_query(parse_blast_tab(blast_raw))
            tool = "blastp_remote_nr"
    else:
        print("No --diamond-db, --blast-db, or --blast-remote-nr; skipping homology search.", file=sys.stderr)

    ann = build_annotation_table(top, seq_map, blast_best, tool)
    ann_path = rep_dir / "stk2_linked_panaroo_annotations.tsv"
    ann.to_csv(ann_path, sep="\t", index=False)
    print(f"Wrote {ann_path}")

    merge_path = rep_dir / "panaroo_display_labels_from_blast.tsv"
    write_label_merge_file(ann, merge_path)
    print(f"Wrote {merge_path} (merge into config/panaroo_display_overrides.tsv manually or cat >>)")

    write_summary_text(
        ann,
        rep_dir / "stk2_linked_panaroo_annotation_summary.txt",
        had_blast=tool not in ("none",),
    )
    print(f"Wrote {rep_dir / 'stk2_linked_panaroo_annotation_summary.txt'}")


if __name__ == "__main__":
    main()
