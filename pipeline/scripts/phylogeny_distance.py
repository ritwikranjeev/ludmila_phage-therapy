#!/usr/bin/env python3
"""Pairwise distances from a Newick tree for pyseer --distances (fixed-effects / MDS).

Based on the pyseer utility (patristic distances via dendropy):
https://github.com/johnlees/pyseer/blob/master/scripts/phylogeny_distance.py

Requires: pip install dendropy
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import List, Set


def _read_phenotype_samples(path: Path) -> List[str]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "sample" not in reader.fieldnames:
            raise SystemExit(
                f"{path}: need a TSV header with a 'sample' column (first column)"
            )
        return [row["sample"].strip() for row in reader if row.get("sample")]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Write square tab-delimited distance matrix (patristic by default)"
    )
    parser.add_argument("phylogeny", help="Tree file (Newick)")
    parser.add_argument("--format", default="newick", help="BioPhylip schema [newick]")
    parser.add_argument(
        "--topology",
        action="store_true",
        help="Ignore branch lengths; use topological path (edge count)",
    )
    parser.add_argument(
        "--phenotypes-tsv",
        type=Path,
        default=None,
        help="Optional: fail if any 'sample' in this file is missing from tree tips",
    )
    args = parser.parse_args()

    try:
        import dendropy
    except ImportError:
        print("ERROR: install dendropy in your pyseer env: pip install dendropy", file=sys.stderr)
        return 1

    tree = dendropy.Tree.get(
        path=args.phylogeny,
        schema=args.format,
        preserve_underscores=True,
    )
    taxa = list(tree.taxon_namespace)
    labels = [t.label for t in taxa]
    tip_set: Set[str] = set(labels)

    if args.phenotypes_tsv is not None:
        needed = _read_phenotype_samples(args.phenotypes_tsv)
        missing = [s for s in needed if s not in tip_set]
        if missing:
            print(
                "ERROR: phenotypes.tsv samples not found as tips in the tree:\n  "
                + "\n  ".join(missing[:40])
                + (f"\n  ... and {len(missing) - 40} more" if len(missing) > 40 else ""),
                file=sys.stderr,
            )
            return 1

    pdm = tree.phylogenetic_distance_matrix()
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow([""] + labels)
    for t1 in taxa:
        row: List[str] = [t1.label]
        for t2 in taxa:
            if args.topology:
                row.append(str(pdm.path_edge_count(t1, t2)))
            else:
                row.append(str(pdm.patristic_distance(t1, t2)))
        w.writerow(row)
    return 0


if __name__ == "__main__":
    sys.exit(main())
