#!/usr/bin/env python3
"""Call defense islands from defense gene table."""

import csv
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, Dict, List, Set, Tuple

PROJECT_ROOT = Path(__file__).resolve().parents[1]
IN_TSV = PROJECT_ROOT / "outputs" / "defense" / "defense_genes.tsv"
OUT_TSV = PROJECT_ROOT / "outputs" / "defense" / "defense_islands.tsv"
OUT_FEATURES = PROJECT_ROOT / "outputs" / "defense" / "defense_family_presence.tsv"
WINDOW_BP = 15000  # editable: cluster genes into islands within this genomic window


def to_int(v: str) -> int:
    try:
        return int(v)
    except Exception:
        return -1


def main() -> None:
    if not IN_TSV.exists():
        raise FileNotFoundError(f"Missing input: {IN_TSV}")
    by_sample_contig: DefaultDict[Tuple[str, str], List[Dict[str, str]]] = defaultdict(list)
    with IN_TSV.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            row["start_i"] = to_int(row.get("start", ""))
            by_sample_contig[(row["sample"], row["contig"])].append(row)

    island_rows: List[Dict[str, str]] = []
    family_presence: DefaultDict[str, Set[str]] = defaultdict(set)
    for (sample, contig), rows in by_sample_contig.items():
        rows = sorted(rows, key=lambda r: r["start_i"])
        island_idx = 0
        cluster = []
        for row in rows:
            family_presence[sample].add(row.get("defense_family", ""))
            if not cluster:
                cluster = [row]
                continue
            if row["start_i"] - cluster[-1]["start_i"] <= WINDOW_BP:
                cluster.append(row)
            else:
                if len(cluster) >= 2:
                    island_idx += 1
                    island_rows.append(make_island(sample, contig, island_idx, cluster))
                cluster = [row]
        if len(cluster) >= 2:
            island_idx += 1
            island_rows.append(make_island(sample, contig, island_idx, cluster))

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_TSV.open("w", newline="") as out:
        fields = [
            "sample", "island_id", "contig", "start", "end", "size_bp",
            "defense_systems", "list_of_locus_tags",
        ]
        w = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(island_rows)

    all_families = sorted({f for fams in family_presence.values() for f in fams if f})
    with OUT_FEATURES.open("w", newline="") as out:
        fields = ["sample"] + all_families
        w = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for sample in sorted(family_presence):
            row = {"sample": sample}
            for fam in all_families:
                row[fam] = 1 if fam in family_presence[sample] else 0
            w.writerow(row)

    print(f"Wrote {OUT_TSV} ({len(island_rows)} islands)")
    print(f"Wrote {OUT_FEATURES}")


def make_island(
    sample: str, contig: str, idx: int, cluster: List[Dict[str, str]]
) -> Dict[str, str]:
    start = min(to_int(r.get("start", "")) for r in cluster)
    end = max(to_int(r.get("end", "")) for r in cluster)
    systems = sorted(
        {r.get("defense_system", "") for r in cluster if r.get("defense_system", "")}
    )
    tags = sorted({r.get("locus_tag", "") for r in cluster if r.get("locus_tag", "")})
    return {
        "sample": sample,
        "island_id": "{}_island_{}".format(sample, idx),
        "contig": contig,
        "start": str(start),
        "end": str(end),
        "size_bp": str(end - start + 1),
        "defense_systems": ",".join(systems),
        "list_of_locus_tags": ",".join(tags),
    }


if __name__ == "__main__":
    main()
