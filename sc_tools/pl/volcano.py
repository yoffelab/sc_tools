"""
Volcano plot utilities.

Provides functions for creating volcano plots with statistical annotations.
"""

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..tl.testing import fdr_correction, mwu


def volcano_plot(
    adata: ad.AnnData,
    group_col: str,
    value_col: str,
    group1: str,
    group2: str,
    ax: plt.Axes | None = None,
    fc_threshold: float = 0.5,
    pval_threshold: float = 0.05,
    annotate_significant: bool = True,
    max_annotations: int | None = None,
) -> tuple[plt.Figure, plt.Axes, pd.DataFrame]:
    """
    Create a volcano plot comparing two groups.

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    group_col : str
        Column name in adata.obs containing group labels
    value_col : str
        Column name in adata.obs containing values to compare
    group1 : str
        First group name
    group2 : str
        Second group name
    ax : Axes, optional
        Matplotlib axes to plot on. If None, creates new figure.
    fc_threshold : float
        Fold change threshold for significance (default: 0.5)
    pval_threshold : float
        P-value threshold for significance (default: 0.05)
    annotate_significant : bool
        Whether to annotate significant points (default: True)
    max_annotations : int, optional
        Maximum number of points to annotate (default: None = all significant)

    Returns
    -------
    tuple
        (figure, axes, results_dataframe)
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    else:
        fig = ax.figure

    # Get values for each group
    group1_values = adata.obs.loc[adata.obs[group_col] == group1, value_col].dropna()
    group2_values = adata.obs.loc[adata.obs[group_col] == group2, value_col].dropna()

    if len(group1_values) < 3 or len(group2_values) < 3:
        ax.text(
            0.5,
            0.5,
            "Insufficient data for comparison",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        return fig, ax, pd.DataFrame()

    # Compute fold change
    mean1 = group1_values.mean()
    mean2 = group2_values.mean()

    if mean1 == 0 and mean2 == 0:
        log2fc = 0.0
    elif mean1 == 0:
        log2fc = np.sign(mean2) * 10
    elif mean2 == 0:
        log2fc = np.sign(mean1) * -10
    else:
        if mean1 > 0 and mean2 > 0:
            log2fc = np.log2(mean2 / mean1)
        else:
            # For z-scores, use difference instead
            log2fc = mean2 - mean1

    # Statistical test
    try:
        stat, pval = mwu(group1_values, group2_values)
    except Exception:
        pval = np.nan

    # Create results dataframe
    results = pd.DataFrame(
        {
            "value_col": [value_col],
            "log2FC": [log2fc],
            "pval": [pval],
            "mean_group1": [mean1],
            "mean_group2": [mean2],
        }
    )

    # Apply FDR correction (if multiple comparisons)
    if pval is not None and not np.isnan(pval):
        _, adj_pvals, _, _ = fdr_correction(np.array([pval]))
        results["adj_pval"] = adj_pvals[0]
    else:
        results["adj_pval"] = np.nan

    # Determine color
    if results["adj_pval"].iloc[0] < pval_threshold and results["log2FC"].iloc[0] > fc_threshold:
        color = "red"
    elif results["adj_pval"].iloc[0] < pval_threshold and results["log2FC"].iloc[0] < -fc_threshold:
        color = "blue"
    elif results["adj_pval"].iloc[0] < pval_threshold:
        color = "orange"
    else:
        color = "gray"

    # Plot
    neg_log10_pval = -np.log10(pval + 1e-300) if not np.isnan(pval) else 0
    ax.scatter(log2fc, neg_log10_pval, s=50, alpha=0.6, color=color)

    # Add threshold lines
    ax.axhline(y=-np.log10(pval_threshold), color="black", linestyle="--", linewidth=1, alpha=0.5)
    ax.axvline(x=fc_threshold, color="black", linestyle="--", linewidth=1, alpha=0.3)
    ax.axvline(x=-fc_threshold, color="black", linestyle="--", linewidth=1, alpha=0.3)

    # Labels
    ax.set_xlabel(f"Log2 Fold Change ({group2} / {group1})", fontsize=11)
    ax.set_ylabel("-Log10(p-value)", fontsize=11)
    ax.set_title(value_col.replace("sig:", "").replace("_z", ""), fontsize=12, fontweight="bold")
    ax.grid(alpha=0.3, linestyle="--")

    return fig, ax, results


def volcano_plot_faceted(
    adata: ad.AnnData,
    group_col: str,
    value_cols: list[str],
    comparisons: list[tuple[str, str, str]],
    figsize: tuple[int, int] = (18, 6),
    fc_threshold: float = 0.5,
    pval_threshold: float = 0.05,
    annotate_significant: bool = True,
    signatures_include: list[str] | None = None,
    signatures_exclude: list[str] | None = None,
) -> tuple[plt.Figure, pd.DataFrame]:
    """
    Create faceted volcano plots for multiple comparisons.

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    group_col : str
        Column name in adata.obs containing group labels
    value_cols : list of str
        List of column names to compare
    comparisons : list of tuples
        List of (group1, group2, comp_name) tuples
    figsize : tuple
        Figure size (default: (18, 6))
    fc_threshold : float
        Fold change threshold (default: 0.5)
    pval_threshold : float
        P-value threshold (default: 0.05)
    annotate_significant : bool
        Whether to annotate significant points (default: True)
    signatures_include : list, optional
        Signatures to include (for filtering)
    signatures_exclude : list, optional
        Signatures to exclude (for filtering)

    Returns
    -------
    tuple
        (figure, results_dataframe)
    """
    # Filter signatures if needed
    if signatures_include is not None or signatures_exclude is not None:
        from ..utils.signatures import filter_signatures

        value_cols = filter_signatures(
            value_cols, include=signatures_include, exclude=signatures_exclude
        )

    # Compute statistics for all comparisons
    all_results = []
    all_plot_data = []

    for group1, group2, comp_name in comparisons:
        results = []

        for value_col in value_cols:
            # Get values
            group1_values = adata.obs.loc[adata.obs[group_col] == group1, value_col].dropna()
            group2_values = adata.obs.loc[adata.obs[group_col] == group2, value_col].dropna()

            if len(group1_values) < 3 or len(group2_values) < 3:
                continue

            # Compute fold change
            mean1 = group1_values.mean()
            mean2 = group2_values.mean()

            if mean1 == 0 and mean2 == 0:
                log2fc = 0.0
            elif mean1 == 0:
                log2fc = np.sign(mean2) * 10
            elif mean2 == 0:
                log2fc = np.sign(mean1) * -10
            else:
                if mean1 > 0 and mean2 > 0:
                    log2fc = np.log2(mean2 / mean1)
                else:
                    log2fc = mean2 - mean1

            # Statistical test
            try:
                stat, pval = mwu(group1_values, group2_values)
            except Exception:
                pval = np.nan

            results.append(
                {
                    "signature": value_col,
                    "log2FC": log2fc,
                    "pval": pval,
                    "mean_group1": mean1,
                    "mean_group2": mean2,
                }
            )

        # Convert to DataFrame and apply FDR
        comp_df = pd.DataFrame(results)
        valid_pvals = comp_df["pval"].dropna()
        if len(valid_pvals) > 0:
            _, adj_pvals, _, _ = fdr_correction(valid_pvals.values)
            comp_df.loc[comp_df["pval"].notna(), "adj_pval"] = adj_pvals
        else:
            comp_df["adj_pval"] = np.nan

        comp_df["comparison"] = comp_name
        comp_df["group1"] = group1
        comp_df["group2"] = group2
        all_results.append(comp_df)

        # Prepare plot data
        plot_df = comp_df.copy()
        plot_df = plot_df.dropna(subset=["log2FC", "pval"])
        plot_df["neg_log10_pval"] = -np.log10(plot_df["pval"] + 1e-300)
        plot_df["color"] = "gray"
        plot_df.loc[
            (plot_df["adj_pval"] < pval_threshold) & (plot_df["log2FC"] > fc_threshold), "color"
        ] = "red"
        plot_df.loc[
            (plot_df["adj_pval"] < pval_threshold) & (plot_df["log2FC"] < -fc_threshold), "color"
        ] = "blue"
        plot_df.loc[
            (plot_df["adj_pval"] < pval_threshold) & (plot_df["log2FC"].abs() <= fc_threshold),
            "color",
        ] = "orange"
        all_plot_data.append(plot_df)

    # Create faceted plot
    fig, axes = plt.subplots(1, len(comparisons), figsize=figsize)
    if len(comparisons) == 1:
        axes = [axes]

    for ax_idx, ((group1, group2, comp_name), plot_df) in enumerate(
        zip(comparisons, all_plot_data, strict=True)
    ):
        ax = axes[ax_idx]

        # Plot points by color
        for color in ["gray", "orange", "red", "blue"]:
            mask = plot_df["color"] == color
            if mask.sum() > 0:
                ax.scatter(
                    plot_df.loc[mask, "log2FC"],
                    plot_df.loc[mask, "neg_log10_pval"],
                    alpha=0.6,
                    s=50,
                    color=color,
                    label={
                        "gray": "Not significant",
                        "orange": f"Significant, |FC| ≤ {fc_threshold}",
                        "red": f"Significant, FC > {fc_threshold}",
                        "blue": f"Significant, FC < -{fc_threshold}",
                    }[color]
                    if ax_idx == 0
                    else None,
                )

        # Add threshold lines
        ax.axhline(
            y=-np.log10(pval_threshold), color="black", linestyle="--", linewidth=1, alpha=0.5
        )
        ax.axvline(x=fc_threshold, color="black", linestyle="--", linewidth=1, alpha=0.3)
        ax.axvline(x=-fc_threshold, color="black", linestyle="--", linewidth=1, alpha=0.3)

        # Labels
        ax.set_xlabel(f"Log2 Fold Change ({group2} / {group1})", fontsize=11)
        ax.set_ylabel("-Log10(p-value)", fontsize=11)
        ax.set_title(comp_name.replace("_", " "), fontsize=12, fontweight="bold")
        ax.grid(alpha=0.3, linestyle="--")

    # Add legend
    if len(comparisons) > 0:
        axes[0].legend(fontsize=8, frameon=True, loc="upper right")

    plt.tight_layout()

    # Combine results
    all_results_df = pd.concat(all_results, ignore_index=True)

    return fig, all_results_df


__all__ = ["volcano_plot", "volcano_plot_faceted"]
