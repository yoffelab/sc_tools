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

    exclude = {"method", "embedding_key", "batch_score", "bio_score", "overall_score"}
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

    exclude = {"method", "embedding_key", "batch_score", "bio_score", "overall_score"}
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

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=df["batch_score"],
            y=df["bio_score"],
            mode="markers+text",
            text=df["method"],
            textposition="top center",
            marker={
                "size": 12,
                "color": df["overall_score"],
                "colorscale": "Viridis",
                "showscale": True,
                "colorbar": {"title": "Overall"},
            },
            hovertemplate=("Method: %{text}<br>Batch: %{x:.3f}<br>Bio: %{y:.3f}<extra></extra>"),
        )
    )
    fig.update_layout(
        title="Batch Removal vs Bio Conservation",
        xaxis_title="Batch Removal Score",
        yaxis_title="Bio Conservation Score",
        xaxis={"range": [-0.05, 1.05]},
        yaxis={"range": [-0.05, 1.05]},
        height=500,
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
