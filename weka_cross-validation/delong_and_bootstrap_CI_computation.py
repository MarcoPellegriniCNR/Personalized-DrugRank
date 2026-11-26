#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd

from sklearn.metrics import roc_auc_score
from scipy.stats import norm


# ============================================================
# DeLong variance for ONE ROC AUC
# ============================================================

def delong_auc_ci(y_true, y_score, confidence=0.95):
    """
    DeLong variance estimate for a single ROC AUC.

    y_true:
        1 = positive class (PD)
        0 = negative class (CR)

    y_score:
        probability/score for PD

    Returns:
        auc
        standard error
        lower CI
        upper CI
    """

    y_true = np.asarray(y_true, dtype=int)
    y_score = np.asarray(y_score, dtype=float)

    pos_scores = y_score[y_true == 1]
    neg_scores = y_score[y_true == 0]

    m = len(pos_scores)
    n = len(neg_scores)

    if m < 2 or n < 2:
        raise ValueError(
            "DeLong CI requires at least two "
            "positive and two negative observations."
        )

    # --------------------------------------------------------
    # Pairwise comparison matrix
    #
    # 1.0 : PD score > CR score
    # 0.5 : tie
    # 0.0 : PD score < CR score
    # --------------------------------------------------------

    comparisons = (
        (pos_scores[:, None] > neg_scores[None, :]).astype(float)
        +
        0.5
        * (pos_scores[:, None] == neg_scores[None, :]).astype(float)
    )

    # Structural components
    v10 = comparisons.mean(axis=1)
    v01 = comparisons.mean(axis=0)

    auc = v10.mean()

    # DeLong variance
    var_auc = (
        np.var(v10, ddof=1) / m
        +
        np.var(v01, ddof=1) / n
    )

    se_auc = np.sqrt(var_auc)

    alpha = 1.0 - confidence

    z = norm.ppf(1.0 - alpha / 2.0)

    lower = auc - z * se_auc
    upper = auc + z * se_auc

    # AUC cannot be outside [0,1]
    lower = max(0.0, lower)
    upper = min(1.0, upper)

    return auc, se_auc, lower, upper


# ============================================================
# Stratified bootstrap
# ============================================================

def stratified_bootstrap_auc_ci(
    y_true,
    y_score,
    n_bootstrap=50000,
    confidence=0.95,
    seed=12345
):
    """
    Bootstrap the SUBJECT-LEVEL held-out scores.

    PD subjects are sampled with replacement from PD subjects.
    CR subjects are sampled with replacement from CR subjects.

    This keeps each bootstrap sample at:
        7 PD + 8 CR
    for the current dataset.
    """

    y_true = np.asarray(y_true, dtype=int)
    y_score = np.asarray(y_score, dtype=float)

    pos_scores = y_score[y_true == 1]
    neg_scores = y_score[y_true == 0]

    n_pos = len(pos_scores)
    n_neg = len(neg_scores)

    rng = np.random.default_rng(seed)

    aucs = np.empty(n_bootstrap, dtype=float)

    for b in range(n_bootstrap):

        boot_pos = rng.choice(
            pos_scores,
            size=n_pos,
            replace=True
        )

        boot_neg = rng.choice(
            neg_scores,
            size=n_neg,
            replace=True
        )

        boot_scores = np.concatenate([
            boot_pos,
            boot_neg
        ])

        boot_labels = np.concatenate([
            np.ones(n_pos, dtype=int),
            np.zeros(n_neg, dtype=int)
        ])

        aucs[b] = roc_auc_score(
            boot_labels,
            boot_scores
        )

    alpha = 1.0 - confidence

    lower = np.quantile(
        aucs,
        alpha / 2.0
    )

    upper = np.quantile(
        aucs,
        1.0 - alpha / 2.0
    )

    return lower, upper, aucs


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        required=True,
        help="CSV containing pooled nested-LOOCV predictions"
    )

    parser.add_argument(
        "--label-col",
        default="binary_actual"
    )

    parser.add_argument(
        "--score-col",
        default="P_PD"
    )

    parser.add_argument(
        "--bootstrap",
        type=int,
        default=50000,
        help="Number of bootstrap replicates"
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=12345
    )

    args = parser.parse_args()

    # --------------------------------------------------------
    # Read data
    # --------------------------------------------------------

    df = pd.read_csv(args.input)

    y_true = df[args.label_col].to_numpy(dtype=int)
    y_score = df[args.score_col].to_numpy(dtype=float)

    n = len(y_true)
    n_pd = np.sum(y_true == 1)
    n_cr = np.sum(y_true == 0)

    if n_pd == 0 or n_cr == 0:
        raise ValueError(
            "Both PD and CR must be present."
        )

    # --------------------------------------------------------
    # Point estimate
    # --------------------------------------------------------

    auc_sklearn = roc_auc_score(
        y_true,
        y_score
    )

    # --------------------------------------------------------
    # DeLong
    # --------------------------------------------------------

    (
        auc_delong,
        se_delong,
        dl_low,
        dl_high
    ) = delong_auc_ci(
        y_true,
        y_score,
        confidence=0.95
    )

    # --------------------------------------------------------
    # Bootstrap
    # --------------------------------------------------------

    (
        boot_low,
        boot_high,
        boot_aucs
    ) = stratified_bootstrap_auc_ci(
        y_true,
        y_score,
        n_bootstrap=args.bootstrap,
        confidence=0.95,
        seed=args.seed
    )

    boot_se = np.std(
        boot_aucs,
        ddof=1
    )

    # --------------------------------------------------------
    # Pair interpretation
    # --------------------------------------------------------

    total_pairs = n_pd * n_cr

    pair_equivalent = auc_sklearn * total_pairs

    # --------------------------------------------------------
    # Output
    # --------------------------------------------------------

    print()
    print("============================================")
    print("NESTED LOOCV POOLED AUC")
    print("============================================")
    print()

    print(f"Subjects             : {n}")
    print(f"PD                    : {n_pd}")
    print(f"CR                    : {n_cr}")
    print(f"PD-CR pairs           : {total_pairs}")
    print()

    print(f"AUC                   : {auc_sklearn:.6f}")
    print(
        f"Equivalent pair credit: "
        f"{pair_equivalent:.3f}/{total_pairs}"
    )

    print()
    print("--------------------------------------------")
    print("DELONG")
    print("--------------------------------------------")

    print(f"AUC                   : {auc_delong:.6f}")
    print(f"DeLong SE             : {se_delong:.6f}")
    print(
        f"95% DeLong CI         : "
        f"[{dl_low:.6f}, {dl_high:.6f}]"
    )

    print()
    print("--------------------------------------------")
    print("STRATIFIED BOOTSTRAP")
    print("--------------------------------------------")

    print(
        f"Bootstrap replicates  : "
        f"{args.bootstrap:,}"
    )

    print(f"Bootstrap SE          : {boot_se:.6f}")

    print(
        f"95% percentile CI     : "
        f"[{boot_low:.6f}, {boot_high:.6f}]"
    )

    print()
    print("============================================")
    print("IMPORTANT")
    print("============================================")
    print(
        "These CIs are conditional on the set of  pooled "
        "outer-LOOCV prediction scores."
    )
    print(
        "They do not fully include variability from "
        "rerunning Auto-WEKA model/feature selection."
    )
    print()


if __name__ == "__main__":
    main()
