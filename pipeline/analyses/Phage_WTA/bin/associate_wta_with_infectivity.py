#!/usr/bin/env python3
"""
Merge WTA feature matrix with Ludmila infectivity and fit a simple logistic regression.

Reads:
  data/wta_features/wta_feature_matrix.tsv
  data/phenotypes/ludmila_infectivity.tsv (lines starting with # are comments)

Writes:
  analysis/wta_infectivity_logreg_coeffs.tsv
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

try:
    import numpy as np
    import pandas as pd
except ImportError as e:
    print(
        "ERROR: numpy and pandas are required. Install in phage_wta_env: pip install numpy pandas",
        file=sys.stderr,
    )
    raise SystemExit(1) from e

try:
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import accuracy_score, roc_auc_score
    from sklearn.model_selection import train_test_split
except ImportError as e:
    print("ERROR: scikit-learn is required. Install in phage_wta_env: pip install scikit-learn", file=sys.stderr)
    raise SystemExit(1) from e

TEST_FRACTION = 0.2
RANDOM_STATE = 42


def project_root() -> Path:
    return Path(__file__).resolve().parent.parent


def load_infectivity(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#", dtype={"host_id": str})


def main() -> int:
    ap = argparse.ArgumentParser(description="WTA features vs infectivity logistic regression.")
    ap.add_argument(
        "--min-hosts",
        type=int,
        default=None,
        metavar="N",
        help="Minimum merged host rows (default: env PHAGE_WTA_ASSOC_MIN_HOSTS or 5). Use 2 for tiny toy tables.",
    )
    args = ap.parse_args()

    min_hosts = args.min_hosts
    if min_hosts is None:
        min_hosts = int(os.environ.get("PHAGE_WTA_ASSOC_MIN_HOSTS", "5"))

    root = project_root()
    wta_path = root / "data" / "wta_features" / "wta_feature_matrix.tsv"
    pheno_path = root / "data" / "phenotypes" / "ludmila_infectivity.tsv"
    out_path = root / "analysis" / "wta_infectivity_logreg_coeffs.tsv"

    if not wta_path.is_file():
        print(f"ERROR: Missing {wta_path}", file=sys.stderr)
        return 1
    if not pheno_path.is_file():
        print(
            f"ERROR: Missing {pheno_path}. Run: python bin/prepare_infectivity_table.py",
            file=sys.stderr,
        )
        return 1

    wta = pd.read_csv(wta_path, sep="\t", dtype={"host_id": str})
    pheno = load_infectivity(pheno_path)

    if "host_id" not in wta.columns or "host_id" not in pheno.columns:
        print("ERROR: host_id column required in both tables.", file=sys.stderr)
        return 1
    if "infectivity" not in pheno.columns:
        print("ERROR: infectivity column required in phenotype table.", file=sys.stderr)
        return 1

    pheno = pheno[["host_id", "infectivity"]].copy()
    pheno["infectivity"] = pd.to_numeric(pheno["infectivity"], errors="coerce")
    pheno = pheno.dropna(subset=["infectivity"])
    pheno = pheno.drop_duplicates(subset=["host_id"], keep="last")

    yr = pheno["infectivity"].to_numpy(dtype=float)
    only01 = np.all(np.isin(yr, (0.0, 1.0)))
    if not only01:
        print(
            "Note: infectivity is not strictly 0/1; binarizing as (infectivity > 0) -> 1.",
            file=sys.stderr,
        )
    y_bin = (yr > 0).astype(int)
    pheno["infectivity_bin"] = y_bin

    merged = wta.merge(pheno[["host_id", "infectivity_bin"]], on="host_id", how="inner")
    if merged.shape[0] < min_hosts:
        print(
            f"ERROR: Too few matched hosts after merge ({merged.shape[0]}; need >= {min_hosts}). "
            "Check host_id spelling, add hosts, or set PHAGE_WTA_ASSOC_MIN_HOSTS / --min-hosts for toy runs.",
            file=sys.stderr,
        )
        return 1

    feature_block = merged.drop(columns=["host_id", "infectivity_bin"])
    X = pd.get_dummies(feature_block, drop_first=False)
    y = merged["infectivity_bin"].values

    if len(np.unique(y)) < 2:
        print(
            "ERROR: Infectivity is constant after merge; need both 0 and 1 (or mixed signs).",
            file=sys.stderr,
        )
        return 1

    strat = y if len(np.unique(y)) > 1 else None
    try:
        X_train, X_test, y_train, y_test = train_test_split(
            X,
            y,
            test_size=TEST_FRACTION,
            random_state=RANDOM_STATE,
            stratify=strat,
        )
    except ValueError:
        X_train, X_test, y_train, y_test = train_test_split(
            X,
            y,
            test_size=TEST_FRACTION,
            random_state=RANDOM_STATE,
            stratify=None,
        )

    clf = LogisticRegression(max_iter=2000, solver="lbfgs")
    clf.fit(X_train, y_train)

    y_pred = clf.predict(X_test)
    acc = accuracy_score(y_test, y_pred)
    try:
        if len(np.unique(y_test)) < 2:
            auc = float("nan")
        else:
            auc = roc_auc_score(y_test, clf.predict_proba(X_test)[:, 1])
    except ValueError:
        auc = float("nan")

    coefs = clf.coef_.ravel()
    names = list(X.columns)
    abs_rank = pd.Series(np.abs(coefs), index=names).rank(ascending=False, method="min").astype(int)

    summary = pd.DataFrame(
        {
            "feature": names,
            "coef": coefs,
            "abs_coef_rank": [abs_rank[n] for n in names],
        }
    ).sort_values("abs_coef_rank")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out_path, sep="\t", index=False)

    print("--- Logistic regression (WTA features -> infectivity) ---")
    print(f"Matched hosts (n): {merged.shape[0]}")
    print(f"Train n={len(y_train)}, test n={len(y_test)}")
    print(f"Test accuracy: {acc:.3f}")
    if np.isnan(auc):
        print("Test ROC AUC: n/a (single class in test split)")
    else:
        print(f"Test ROC AUC: {auc:.3f}")
    print("\nCoefficients (by |coef| rank):")
    for _, row in summary.sort_values("abs_coef_rank").iterrows():
        print(f"  {row['abs_coef_rank']:3d}  {row['feature']:40s}  {row['coef']:+.4f}")
    print(f"\nWrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
