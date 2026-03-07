"""
Benchmark plotting functions for segmentation and integration comparison.

Functions return ``matplotlib.Figure`` for static plots (UMAP, overlay)
or ``plotly.graph_objects.Figure`` for interactive charts (used in HTML reports).
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from anndata import AnnData

__all__ = [
    # Segmentation
    "plot_segmentation_comparison_table",
    "plot_morphology_distributions",
    "plot_marker_snr_heatmap",
    "plot_segmentation_overlay",
    "plot_composite_score_bar",
    # Integration
    "plot_integration_comparison_table",
    "plot_integration_radar",
    "plot_integration_ranking_bar",
    "plot_batch_vs_bio",
    "plot_embedding_comparison_umap",
    # Extended benchmark
    "plot_strategy_comparison",
    "plot_dataset_heatmap",
    "plot_tissue_boxplot",
    "plot_generalization_radar",
    "plot_runtime_comparison",
    "plot_failure_gallery",
    "plot_metric_correlation_scatter",
]

logger = logging.getLogger(__name__)


def _require_plotly():
    """Import and return plotly.graph_objects, raising helpful error if missing."""
    try:
        import plotly.graph_objects as go

        return go
    except ImportError as e:
        raise ImportError(
            "plotly is required for interactive benchmark plots. "
            "Install with: pip install sc-tools[benchmark]"
        ) from e


# ---------------------------------------------------------------------------
# Segmentation plots
# ---------------------------------------------------------------------------


def plot_segmentation_comparison_table(df: pd.DataFrame):
    """Heatmap table of segmentation comparison, best values highlighted.

    Parameters
    ----------
    df
        Output from ``compare_segmentations``.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    go = _require_plotly()

    display_cols = [c for c in df.columns if c not in ("method",)]
    methods = df["method"].tolist()
    z = df[display_cols].values.astype(float)

    fig = go.Figure(
        data=go.Heatmap(
            z=z,
            x=display_cols,
            y=methods,
            colorscale="RdYlGn",
            text=np.round(z, 2).astype(str),
            texttemplate="%{text}",
            hovertemplate="Method: %{y}<br>Metric: %{x}<br>Value: %{z:.3f}<extra></extra>",
        )
    )
    fig.update_layout(
        title="Segmentation Comparison",
        xaxis_title="Metric",
        yaxis_title="Method",
        height=max(300, 60 * len(methods)),
    )
    return fig


def plot_morphology_distributions(
    morphology_dict: dict[str, pd.DataFrame],
    metrics: list[str] | None = None,
):
    """Overlapping violin/KDE per metric per method.

    Parameters
    ----------
    morphology_dict
        Dict mapping method name to morphology DataFrame
        (output of ``compute_morphology_metrics``).
    metrics
        Which columns to plot. Default: circularity, solidity, eccentricity.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    go = _require_plotly()

    if metrics is None:
        metrics = ["circularity", "solidity", "eccentricity"]

    from plotly.subplots import make_subplots

    fig = make_subplots(rows=1, cols=len(metrics), subplot_titles=metrics)

    for col_idx, metric in enumerate(metrics, 1):
        for method_name, morph_df in morphology_dict.items():
            if metric in morph_df.columns:
                fig.add_trace(
                    go.Violin(
                        y=morph_df[metric].dropna(),
                        name=method_name,
                        legendgroup=method_name,
                        showlegend=(col_idx == 1),
                        scalemode="width",
                        width=0.8,
                    ),
                    row=1,
                    col=col_idx,
                )

    fig.update_layout(
        title="Morphology Distributions by Method",
        height=400,
        violinmode="overlay",
    )
    return fig


def plot_marker_snr_heatmap(marker_quality_dict: dict[str, pd.DataFrame]):
    """Heatmap: rows=markers, cols=methods, values=SNR.

    Parameters
    ----------
    marker_quality_dict
        Dict mapping method name to marker quality DataFrame
        (output of ``compute_marker_quality``).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    go = _require_plotly()

    methods = list(marker_quality_dict.keys())
    all_markers = set()
    for mq in marker_quality_dict.values():
        all_markers.update(mq["marker"].tolist())
    all_markers = sorted(all_markers)

    z = np.zeros((len(all_markers), len(methods)))
    for j, method in enumerate(methods):
        mq = marker_quality_dict[method].set_index("marker")
        for i, marker in enumerate(all_markers):
            if marker in mq.index:
                z[i, j] = mq.loc[marker, "snr"]

    fig = go.Figure(
        data=go.Heatmap(
            z=z,
            x=methods,
            y=all_markers,
            colorscale="Viridis",
            text=np.round(z, 2).astype(str),
            texttemplate="%{text}",
            hovertemplate="Marker: %{y}<br>Method: %{x}<br>SNR: %{z:.2f}<extra></extra>",
        )
    )
    fig.update_layout(
        title="Marker SNR by Method",
        xaxis_title="Method",
        yaxis_title="Marker",
        height=max(300, 30 * len(all_markers)),
    )
    return fig


def plot_segmentation_overlay(
    intensity: np.ndarray,
    masks: dict[str, np.ndarray],
    roi: tuple[int, int, int, int] | None = None,
) -> plt.Figure:
    """Side-by-side boundary overlays on intensity image.

    Parameters
    ----------
    intensity
        2D intensity image (grayscale).
    masks
        Dict mapping method name to labeled mask.
    roi
        Optional (row_start, row_end, col_start, col_end) crop.

    Returns
    -------
    matplotlib.Figure
    """
    from skimage.segmentation import find_boundaries

    n_methods = len(masks)
    fig, axes = plt.subplots(1, n_methods, figsize=(5 * n_methods, 5))
    if n_methods == 1:
        axes = [axes]

    for ax, (name, mask) in zip(axes, masks.items(), strict=False):
        img = intensity.copy()
        m = mask.copy()
        if roi is not None:
            r0, r1, c0, c1 = roi
            img = img[r0:r1, c0:c1]
            m = m[r0:r1, c0:c1]

        ax.imshow(img, cmap="gray")
        boundaries = find_boundaries(m, mode="thick")
        overlay = np.zeros((*img.shape[:2], 4))
        overlay[boundaries] = [1, 0, 0, 0.7]  # red boundaries
        ax.imshow(overlay)
        ax.set_title(name)
        ax.axis("off")

    fig.tight_layout()
    return fig


def plot_composite_score_bar(df: pd.DataFrame):
    """Horizontal bar chart ranked by composite score.

    Parameters
    ----------
    df
        Output from ``compare_segmentations``.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    go = _require_plotly()

    sorted_df = df.sort_values("composite_score", ascending=True)

    fig = go.Figure(
        go.Bar(
            x=sorted_df["composite_score"],
            y=sorted_df["method"],
            orientation="h",
            marker_color="steelblue",
            text=sorted_df["composite_score"].round(1),
            textposition="auto",
            hovertemplate="Method: %{y}<br>Score: %{x:.1f}<extra></extra>",
        )
    )
    fig.update_layout(
        title="Composite Segmentation Score",
        xaxis_title="Score (0-100)",
        yaxis_title="Method",
        height=max(300, 50 * len(df)),
        xaxis={"range": [0, 105]},
    )
    return fig


# ---------------------------------------------------------------------------
# Integration plots
# ---------------------------------------------------------------------------


def plot_integration_comparison_table(df: pd.DataFrame):
    """Heatmap table with batch/bio metric headers.

    Parameters
    ----------
    df
        Output from ``compare_integrations``.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    go = _require_plotly()

    exclude = {"method", "embedding_key", "batch_score", "bio_score", "overall_score", "group"}
    metric_cols = [c for c in df.columns if c not in exclude]
    methods = df["method"].tolist()
    z = df[metric_cols].values.astype(float)

    fig = go.Figure(
        data=go.Heatmap(
            z=z,
            x=metric_cols,
            y=methods,
            colorscale="RdYlGn",
            zmin=0,
            zmax=1,
            text=np.round(z, 3).astype(str),
            texttemplate="%{text}",
            hovertemplate="Method: %{y}<br>Metric: %{x}<br>Value: %{z:.3f}<extra></extra>",
        )
    )
    fig.update_layout(
        title="Integration Metrics Comparison",
        xaxis_title="Metric",
        yaxis_title="Method",
        height=max(300, 60 * len(methods)),
    )
    return fig


def plot_integration_radar(df: pd.DataFrame):
    """Radar/spider plot per integration method.

    Parameters
    ----------
    df
        Output from ``compare_integrations``.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    go = _require_plotly()

    exclude = {"method", "embedding_key", "batch_score", "bio_score", "overall_score", "group"}
    metric_cols = [c for c in df.columns if c not in exclude]

    fig = go.Figure()
    for _, row in df.iterrows():
        values = [row[c] for c in metric_cols]
        values.append(values[0])  # close the polygon
        categories = metric_cols + [metric_cols[0]]

        fig.add_trace(
            go.Scatterpolar(
                r=values,
                theta=categories,
                name=row["method"],
                fill="toself",
                opacity=0.6,
            )
        )

    fig.update_layout(
        title="Integration Metrics Radar",
        polar={"radialaxis": {"visible": True, "range": [0, 1]}},
        height=500,
    )
    return fig


def plot_integration_ranking_bar(df: pd.DataFrame):
    """Grouped bar: overall/batch/bio score per method.

    Parameters
    ----------
    df
        Output from ``compare_integrations``.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    go = _require_plotly()

    sorted_df = df.sort_values("overall_score", ascending=False)

    fig = go.Figure()
    for score_col, color, label in [
        ("overall_score", "#2c3e50", "Overall"),
        ("batch_score", "#e74c3c", "Batch Removal"),
        ("bio_score", "#27ae60", "Bio Conservation"),
    ]:
        if score_col in sorted_df.columns:
            fig.add_trace(
                go.Bar(
                    x=sorted_df["method"],
                    y=sorted_df[score_col],
                    name=label,
                    marker_color=color,
                    text=sorted_df[score_col].round(3),
                    textposition="auto",
                )
            )

    fig.update_layout(
        title="Integration Ranking",
        xaxis_title="Method",
        yaxis_title="Score",
        barmode="group",
        yaxis={"range": [0, 1.05]},
        height=400,
    )
    return fig


def plot_batch_vs_bio(df: pd.DataFrame):
    """Scatter: x=batch_score, y=bio_score, labeled points.

    Parameters
    ----------
    df
        Output from ``compare_integrations``.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    go = _require_plotly()

    # Spread text labels to reduce overlap
    _positions = [
        "top center",
        "bottom center",
        "top right",
        "bottom left",
        "top left",
        "bottom right",
        "middle right",
        "middle left",
    ]
    n = len(df)
    text_pos = [_positions[i % len(_positions)] for i in range(n)]

    # Convert to plain Python lists to avoid binary bdata serialization
    # (bdata not supported by plotly.js 2.27 CDN)
    x_vals = df["batch_score"].tolist()
    y_vals = df["bio_score"].tolist()
    color_vals = df["overall_score"].tolist()
    labels = df["method"].tolist()

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=x_vals,
            y=y_vals,
            mode="markers+text",
            text=labels,
            textposition=text_pos,
            marker={
                "size": 14,
                "color": color_vals,
                "colorscale": "Viridis",
                "showscale": True,
                "colorbar": {"title": "Overall"},
            },
            hovertemplate=("Method: %{text}<br>Batch: %{x:.3f}<br>Bio: %{y:.3f}<extra></extra>"),
            textfont={"size": 11},
        )
    )
    fig.update_layout(
        title="Batch Removal vs Bio Conservation",
        xaxis_title="Batch Removal Score",
        yaxis_title="Bio Conservation Score",
        xaxis={"range": [-0.05, 1.05]},
        yaxis={"range": [-0.05, 1.05]},
        height=550,
    )
    return fig


def plot_embedding_comparison_umap(
    adata: AnnData,
    embeddings: dict[str, str],
    color_by: str,
    n_cols: int = 3,
) -> plt.Figure:
    """Grid of UMAPs colored by a given obs column.

    Parameters
    ----------
    adata
        AnnData with embeddings in ``obsm``.
    embeddings
        Dict mapping method name to ``obsm`` key.
    color_by
        Column in ``adata.obs`` to color by.
    n_cols
        Number of columns in the grid.

    Returns
    -------
    matplotlib.Figure
    """
    import scanpy as sc

    n_methods = len(embeddings)
    n_rows = max(1, (n_methods + n_cols - 1) // n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))

    if n_methods == 1:
        axes = np.array([axes])
    axes = np.atleast_2d(axes)

    for idx, (name, key) in enumerate(embeddings.items()):
        row, col = divmod(idx, n_cols)
        ax = axes[row, col]

        if key not in adata.obsm:
            ax.set_title(f"{name} (missing)")
            ax.axis("off")
            continue

        # Compute UMAP for this embedding
        tmp = adata.copy()
        sc.pp.neighbors(tmp, use_rep=key, n_neighbors=15)
        sc.tl.umap(tmp)

        sc.pl.umap(tmp, color=color_by, ax=ax, show=False, title=name, frameon=False)

    # Turn off unused axes
    for idx in range(n_methods, n_rows * n_cols):
        row, col = divmod(idx, n_cols)
        axes[row, col].axis("off")

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Extended benchmark plots
# ---------------------------------------------------------------------------


def plot_strategy_comparison(
    results_df: pd.DataFrame,
    metric: str = "n_cells",
):
    """Box plot comparing strategies across all ROIs.

    Parameters
    ----------
    results_df
        Full benchmark results with ``strategy`` and metric columns.
    metric
        Metric to compare.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    go = _require_plotly()

    fig = go.Figure()
    for strategy in sorted(results_df["strategy"].unique()):
        subset = results_df[results_df["strategy"] == strategy]
        if metric in subset.columns:
            fig.add_trace(
                go.Box(
                    y=subset[metric],
                    name=f"Strategy {strategy}",
                    boxpoints="outliers",
                )
            )

    fig.update_layout(
        title=f"Strategy Comparison: {metric}",
        yaxis_title=metric,
        height=400,
    )
    return fig


def plot_dataset_heatmap(
    results_df: pd.DataFrame,
    metric: str = "n_cells",
):
    """Heatmap: methods x datasets, values = mean metric.

    Parameters
    ----------
    results_df
        Full benchmark results.
    metric
        Metric to display.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    go = _require_plotly()

    if metric not in results_df.columns:
        fig = go.Figure()
        fig.update_layout(title=f"Dataset Heatmap: {metric} (no data)")
        return fig

    pivot = results_df.pivot_table(
        values=metric,
        index="method",
        columns="dataset",
        aggfunc="mean",
    )

    fig = go.Figure(
        data=go.Heatmap(
            z=pivot.values,
            x=pivot.columns.tolist(),
            y=pivot.index.tolist(),
            colorscale="RdYlGn",
            text=np.round(pivot.values, 1).astype(str),
            texttemplate="%{text}",
            hovertemplate="Method: %{y}<br>Dataset: %{x}<br>Value: %{z:.2f}<extra></extra>",
        )
    )
    fig.update_layout(
        title=f"Per-Dataset Performance: {metric}",
        xaxis_title="Dataset",
        yaxis_title="Method",
        height=max(300, 40 * len(pivot)),
    )
    return fig


def plot_tissue_boxplot(
    results_df: pd.DataFrame,
    metric: str = "n_cells",
):
    """Box plots per tissue type, colored by method.

    Parameters
    ----------
    results_df
        Full benchmark results.
    metric
        Metric to display.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    go = _require_plotly()

    if "tissue" not in results_df.columns or metric not in results_df.columns:
        fig = go.Figure()
        fig.update_layout(title=f"Tissue Analysis: {metric} (no data)")
        return fig

    fig = go.Figure()
    for method in sorted(results_df["method"].unique()):
        subset = results_df[results_df["method"] == method]
        fig.add_trace(
            go.Box(
                x=subset["tissue"],
                y=subset[metric],
                name=method,
                boxpoints=False,
            )
        )

    fig.update_layout(
        title=f"Tissue-Specific Performance: {metric}",
        xaxis_title="Tissue",
        yaxis_title=metric,
        boxmode="group",
        height=400,
    )
    return fig


def plot_generalization_radar(
    gen_matrix: pd.DataFrame,
    top_n: int = 5,
):
    """Radar plot showing method generalization across datasets.

    Parameters
    ----------
    gen_matrix
        Output from ``compute_generalization_matrix()``.
    top_n
        Number of top methods to show.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    go = _require_plotly()

    if "overall_mean" in gen_matrix.columns:
        top = gen_matrix.nlargest(top_n, "overall_mean")
        datasets = [c for c in gen_matrix.columns if c != "overall_mean"]
    else:
        top = gen_matrix.head(top_n)
        datasets = gen_matrix.columns.tolist()

    fig = go.Figure()
    for method, row in top.iterrows():
        values = [row[d] for d in datasets]
        values.append(values[0])  # close polygon
        cats = datasets + [datasets[0]]

        fig.add_trace(
            go.Scatterpolar(
                r=values,
                theta=cats,
                name=str(method),
                fill="toself",
                opacity=0.6,
            )
        )

    fig.update_layout(
        title="Generalization Across Datasets",
        polar={"radialaxis": {"visible": True}},
        height=500,
    )
    return fig


def plot_runtime_comparison(
    results_df: pd.DataFrame,
):
    """Bar chart of mean runtime per method.

    Parameters
    ----------
    results_df
        Full benchmark results with ``runtime_s`` column.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    go = _require_plotly()

    if "runtime_s" not in results_df.columns:
        fig = go.Figure()
        fig.update_layout(title="Runtime Comparison (no data)")
        return fig

    runtime = results_df.groupby("method")["runtime_s"].agg(["mean", "std"]).reset_index()
    runtime = runtime.sort_values("mean")

    fig = go.Figure(
        go.Bar(
            x=runtime["mean"],
            y=runtime["method"],
            orientation="h",
            marker_color="steelblue",
            error_x={"type": "data", "array": runtime["std"], "visible": True},
            text=runtime["mean"].round(1),
            textposition="auto",
        )
    )
    fig.update_layout(
        title="Runtime Comparison (seconds per ROI)",
        xaxis_title="Runtime (s)",
        yaxis_title="Method",
        height=max(300, 40 * len(runtime)),
    )
    return fig


def plot_failure_gallery(
    results_df: pd.DataFrame,
    intensity_images: dict[str, np.ndarray] | None = None,
    masks: dict[str, dict[str, np.ndarray]] | None = None,
    n_worst: int = 3,
    metric: str = "boundary_regularity",
) -> plt.Figure:
    """Show worst ROIs per method as an image grid.

    Parameters
    ----------
    results_df
        Full benchmark results.
    intensity_images
        Dict of roi_id -> intensity image.
    masks
        Dict of roi_id -> {method: mask}.
    n_worst
        Number of worst ROIs per method to show.
    metric
        Metric to rank by (ascending = worst first).

    Returns
    -------
    matplotlib.Figure (or empty figure if no images available).
    """
    if intensity_images is None or masks is None:
        fig, ax = plt.subplots(1, 1, figsize=(6, 2))
        ax.text(
            0.5,
            0.5,
            "No image data available for failure gallery",
            ha="center",
            va="center",
            fontsize=12,
        )
        ax.axis("off")
        return fig

    from sc_tools.bm.analysis import identify_failure_cases

    worst = identify_failure_cases(results_df, metric=metric, n_worst=n_worst)
    methods = worst["method"].unique()

    n_m = len(methods)
    fig, axes = plt.subplots(n_m, n_worst, figsize=(4 * n_worst, 4 * n_m))
    if n_m == 1:
        axes = axes[np.newaxis, :]

    for i, method in enumerate(methods):
        method_worst = worst[worst["method"] == method].head(n_worst)
        for j, (_, row) in enumerate(method_worst.iterrows()):
            ax = axes[i, j]
            roi_id = row.get("roi_id", "")
            if roi_id in intensity_images:
                ax.imshow(intensity_images[roi_id], cmap="gray")
            ax.set_title(f"{method}\n{roi_id}", fontsize=8)
            ax.axis("off")

    fig.tight_layout()
    return fig


def plot_metric_correlation_scatter(
    results_df: pd.DataFrame,
    metric_x: str = "n_cells",
    metric_y: str = "boundary_regularity",
):
    """Scatter plot of two metrics colored by method.

    Parameters
    ----------
    results_df
        Full benchmark results.
    metric_x
        X-axis metric.
    metric_y
        Y-axis metric.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    go = _require_plotly()

    if metric_x not in results_df.columns or metric_y not in results_df.columns:
        fig = go.Figure()
        fig.update_layout(title=f"{metric_x} vs {metric_y} (no data)")
        return fig

    fig = go.Figure()
    for method in sorted(results_df["method"].unique()):
        subset = results_df[results_df["method"] == method]
        fig.add_trace(
            go.Scatter(
                x=subset[metric_x],
                y=subset[metric_y],
                mode="markers",
                name=method,
                opacity=0.6,
                hovertemplate=(
                    f"Method: {method}<br>"
                    f"{metric_x}: %{{x:.2f}}<br>"
                    f"{metric_y}: %{{y:.2f}}"
                    "<extra></extra>"
                ),
            )
        )

    fig.update_layout(
        title=f"{metric_x} vs {metric_y}",
        xaxis_title=metric_x,
        yaxis_title=metric_y,
        height=500,
    )
    return fig
