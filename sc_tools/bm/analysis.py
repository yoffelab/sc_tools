"""Cross-dataset statistics and generalization analysis for benchmarks.

Computes generalization matrices, tissue-specific rankings, statistical
comparisons (paired Wilcoxon + BH correction), and identifies failure cases.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from scipy import stats

__all__ = [
    "compute_generalization_matrix",
    "rank_methods_per_tissue",
    "statistical_comparison",
    "identify_failure_cases",
]

logger = logging.getLogger(__name__)


def compute_generalization_matrix(
    results_df: pd.DataFrame,
    metric: str = "boundary_regularity",
) -> pd.DataFrame:
    """Compute cross-dataset performance matrix per method.

    Parameters
    ----------
    results_df
        Full benchmark results.
    metric
        Metric column to summarize.

    Returns
    -------
    DataFrame with rows=methods, columns=datasets, values=mean metric.
    """
    if metric not in results_df.columns:
        available = [c for c in results_df.columns if results_df[c].dtype in (np.float64, np.int64)]
        raise ValueError(f"Metric {metric!r} not found. Available: {available}")

    pivot = results_df.pivot_table(
        values=metric,
        index="method",
        columns="dataset",
        aggfunc="mean",
    )

    # Add overall mean column
    pivot["overall_mean"] = pivot.mean(axis=1)
    pivot = pivot.sort_values("overall_mean", ascending=False)

    return pivot


def rank_methods_per_tissue(
    results_df: pd.DataFrame,
    metric: str = "boundary_regularity",
) -> pd.DataFrame:
    """Rank methods within each tissue type.

    Parameters
    ----------
    results_df
        Full benchmark results.
    metric
        Metric to rank by.

    Returns
    -------
    DataFrame with columns: tissue, method, mean_{metric}, rank.
    """
    if "tissue" not in results_df.columns or metric not in results_df.columns:
        return pd.DataFrame()

    grouped = (
        results_df.groupby(["tissue", "method"])[metric].agg(["mean", "std", "count"]).reset_index()
    )
    grouped.columns = ["tissue", "method", f"mean_{metric}", f"std_{metric}", "n_rois"]

    # Rank within each tissue
    grouped["rank"] = grouped.groupby("tissue")[f"mean_{metric}"].rank(
        ascending=False, method="min"
    )
    grouped = grouped.sort_values(["tissue", "rank"])

    return grouped


def statistical_comparison(
    results_df: pd.DataFrame,
    metric: str = "boundary_regularity",
    group_col: str = "method",
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Pairwise Wilcoxon signed-rank tests with BH correction.

    Parameters
    ----------
    results_df
        Full benchmark results (must have per-ROI values).
    metric
        Metric to compare.
    group_col
        Column defining groups (methods).
    alpha
        Significance threshold.

    Returns
    -------
    DataFrame with columns: method_a, method_b, statistic, p_value,
    p_adjusted, significant.
    """
    if metric not in results_df.columns:
        return pd.DataFrame()

    # Create paired data: need common ROIs between methods
    methods = results_df[group_col].unique()
    roi_col = "roi_id" if "roi_id" in results_df.columns else None

    if roi_col is None:
        logger.warning("No roi_id column; cannot perform paired tests")
        return pd.DataFrame()

    comparisons = []
    for i, m_a in enumerate(methods):
        for m_b in methods[i + 1 :]:
            df_a = results_df[results_df[group_col] == m_a].set_index(roi_col)
            df_b = results_df[results_df[group_col] == m_b].set_index(roi_col)

            # Find common ROIs
            common = df_a.index.intersection(df_b.index)
            if len(common) < 5:
                continue

            vals_a = df_a.loc[common, metric].values
            vals_b = df_b.loc[common, metric].values

            try:
                stat, p_val = stats.wilcoxon(vals_a, vals_b, alternative="two-sided")
            except ValueError:
                continue

            comparisons.append(
                {
                    "method_a": m_a,
                    "method_b": m_b,
                    "statistic": float(stat),
                    "p_value": float(p_val),
                    "n_pairs": len(common),
                    "mean_diff": float(np.mean(vals_a - vals_b)),
                }
            )

    if not comparisons:
        return pd.DataFrame()

    df = pd.DataFrame(comparisons)

    # Benjamini-Hochberg correction
    n = len(df)
    df = df.sort_values("p_value")
    df["rank"] = range(1, n + 1)
    df["p_adjusted"] = df["p_value"] * n / df["rank"]
    df["p_adjusted"] = df["p_adjusted"].clip(upper=1.0)

    # Ensure monotonicity
    df["p_adjusted"] = df["p_adjusted"][::-1].cummin()[::-1]

    df["significant"] = df["p_adjusted"] < alpha
    df = df.drop(columns=["rank"]).sort_values("p_adjusted")

    return df


def identify_failure_cases(
    results_df: pd.DataFrame,
    metric: str = "boundary_regularity",
    n_worst: int = 10,
) -> pd.DataFrame:
    """Identify worst-performing ROIs per method.

    Parameters
    ----------
    results_df
        Full benchmark results.
    metric
        Metric to use (lower = worse).
    n_worst
        Number of worst ROIs to return per method.

    Returns
    -------
    DataFrame with worst ROIs per method.
    """
    if metric not in results_df.columns:
        return pd.DataFrame()

    worst = (
        results_df.sort_values(metric, ascending=True)
        .groupby("method")
        .head(n_worst)
        .sort_values(["method", metric])
    )

    return worst
