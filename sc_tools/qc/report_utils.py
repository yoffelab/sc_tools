"""
Shared utilities for QC HTML report generation.

Provides figure encoding, template rendering, date stamping, metric table
builders, embedding auto-detection, and optional segmentation/integration
section computation.
"""

from __future__ import annotations

import base64
import io
import logging
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData

__all__ = [
    "fig_to_base64",
    "plotly_to_html",
    "render_template",
    "get_date_stamp",
    "get_modality_terms",
    "build_metrics_table_rows",
    "auto_detect_embeddings",
    "compute_segmentation_section",
    "compute_integration_section",
]

logger = logging.getLogger(__name__)

# Modalities that use protein markers instead of gene expression
_PROTEIN_MODALITIES = {"imc"}

# Modalities that produce single-cell (not spot) observations
_CELL_MODALITIES = {
    "imc",
    "xenium",
    "cosmx",
    "cosmx_1k",
    "cosmx_6k",
    "cosmx_full_library",
    "visium_hd_cell",
}


def get_modality_terms(modality: str) -> dict[str, str | bool]:
    """Return display-label terminology adapted to the data modality.

    Internal column names (``n_genes_by_counts``, ``total_counts``) are
    unchanged -- only **display labels** in plots and HTML templates differ.

    Parameters
    ----------
    modality
        One of ``visium``, ``visium_hd``, ``xenium``, ``cosmx``, ``imc``, etc.

    Returns
    -------
    dict
        Keys: ``feature``, ``features``, ``feature_lower``, ``features_lower``,
        ``observation``, ``observations``, ``observation_lower``,
        ``observations_lower``, ``intensity_label``, ``feature_count_label``,
        ``has_mt``.
    """
    is_protein = modality in _PROTEIN_MODALITIES
    is_cell = modality in _CELL_MODALITIES
    return {
        "feature": "Protein" if is_protein else "Gene",
        "features": "Proteins" if is_protein else "Genes",
        "feature_lower": "protein" if is_protein else "gene",
        "features_lower": "proteins" if is_protein else "genes",
        "observation": "Cell" if is_cell else "Spot",
        "observations": "Cells" if is_cell else "Spots",
        "observation_lower": "cell" if is_cell else "spot",
        "observations_lower": "cells" if is_cell else "spots",
        "intensity_label": "Total intensity" if is_protein else "Total counts",
        "feature_count_label": "Proteins detected" if is_protein else "Genes detected",
        "has_mt": not is_protein,
    }


# Known integration embedding keys in priority order
_KNOWN_EMBEDDINGS = {
    "scVI": "X_scVI",
    "scVI (raw)": "X_scVI_raw",
    "scANVI": "X_scANVI",
    "scANVI (raw)": "X_scANVI_raw",
    "Harmony": "X_pca_harmony",
    "CytoVI": "X_cytovi",
    "CytoVI (arcsinh)": "X_cytovi",
    "BBKNN": "X_umap_bbknn",
    "ComBat": "X_pca_combat",
    "Scanorama": "X_scanorama",
    "DestVI": "X_destvi",
    "IMC Phenotyping": "X_pca_imc_phenotyping",
    "IMC Pheno + Harmony": "X_imc_pheno_harmony",
    "Z-score + Harmony": "X_zscore_harmony",
    "Unintegrated (PCA)": "X_pca",
}


def fig_to_base64(fig: plt.Figure, dpi: int = 150) -> str:
    """Render a matplotlib figure to a base64-encoded PNG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("ascii")


def plotly_to_html(fig) -> str:
    """Convert a plotly figure to an embeddable HTML div string.

    Parameters
    ----------
    fig
        A ``plotly.graph_objects.Figure``.

    Returns
    -------
    str
        HTML ``<div>`` with embedded plotly chart (no full page wrapper).
    """
    return fig.to_html(full_html=False, include_plotlyjs=False)


def render_template(template_path: str | Path, context: dict) -> str:
    """Render a Jinja2 template with the given context.

    Parameters
    ----------
    template_path
        Path to the ``.html`` template file.
    context
        Dict of variables passed to ``Template.render()``.

    Returns
    -------
    str
        Rendered HTML string.
    """
    try:
        from jinja2 import Template
    except ImportError as e:
        raise ImportError(
            "jinja2 is required for HTML report generation. Install with: pip install jinja2"
        ) from e

    template_text = Path(template_path).read_text()
    template = Template(template_text)
    return template.render(**context)


def get_date_stamp(date_stamp: str | None = None) -> str:
    """Return a YYYYMMDD date stamp (today if *date_stamp* is ``None``)."""
    if date_stamp is not None:
        return date_stamp
    return datetime.now().strftime("%Y%m%d")


def build_metrics_table_rows(classified: pd.DataFrame) -> list[dict]:
    """Build per-sample table rows from a ``classify_samples`` DataFrame.

    Parameters
    ----------
    classified
        Output of ``sc_tools.qc.classify_samples``.

    Returns
    -------
    list[dict]
        Each dict has keys: ``sample``, ``qc_pass``, ``qc_fail_reasons``,
        ``n_spots``, ``total_counts_median``, ``n_genes_median``,
        ``n_genes_detected``, ``pct_mt_median``, ``pct_mt_gt5``.
    """
    table_rows: list[dict] = []
    for sample in classified.index:
        row: dict = {"sample": str(sample)}
        row["qc_pass"] = bool(classified.loc[sample, "qc_pass"])
        row["qc_fail_reasons"] = (
            str(classified.loc[sample, "qc_fail_reasons"])
            if classified.loc[sample, "qc_fail_reasons"]
            else ""
        )
        for col in [
            "n_spots",
            "total_counts_median",
            "n_genes_median",
            "n_genes_detected",
            "pct_mt_median",
            "pct_mt_gt5",
        ]:
            val = classified.loc[sample, col] if col in classified.columns else None
            if val is not None and not (isinstance(val, float) and np.isnan(val)):
                row[col] = val
            else:
                row[col] = None
        table_rows.append(row)
    return table_rows


def auto_detect_embeddings(adata: AnnData) -> dict[str, str]:
    """Scan ``adata.obsm`` for known integration embedding keys.

    Returns
    -------
    dict[str, str]
        Mapping of human-readable name to ``obsm`` key, for keys that exist.
    """
    found: dict[str, str] = {}
    for name, key in _KNOWN_EMBEDDINGS.items():
        if key in adata.obsm:
            found[name] = key
    return found


def compute_segmentation_section(
    adata: AnnData,
    masks_dir: str | Path,
    sample_col: str = "library_id",
) -> dict | None:
    """Load masks and compute per-ROI segmentation scores.

    Parameters
    ----------
    adata
        AnnData with sample annotations.
    masks_dir
        Directory containing mask TIFF files (one per ROI/sample).
    sample_col
        Column in ``adata.obs`` identifying samples.

    Returns
    -------
    dict or None
        Dict with ``"df"`` (DataFrame of per-ROI scores) and ``"plots"``
        (dict of plot name to base64 string). Returns ``None`` if no masks
        are found.
    """
    masks_dir = Path(masks_dir)
    if not masks_dir.exists():
        logger.info("Segmentation masks dir not found: %s", masks_dir)
        return None

    mask_files = sorted(masks_dir.glob("*.tif")) + sorted(masks_dir.glob("*.tiff"))
    if not mask_files:
        logger.info("No mask TIFFs found in %s", masks_dir)
        return None

    try:
        from sc_tools.bm.mask_io import load_mask
        from sc_tools.bm.segmentation import score_segmentation
    except ImportError:
        logger.warning("sc_tools.bm not available; skipping segmentation scoring")
        return None

    rows = []
    for mf in mask_files:
        try:
            mask = load_mask(str(mf))
            scores = score_segmentation(mask)
            scores["roi"] = mf.stem
            rows.append(scores)
        except Exception:
            logger.warning("Failed to score mask %s", mf.name, exc_info=True)
    if not rows:
        return None

    df = pd.DataFrame(rows).set_index("roi")

    # Composite score bar chart
    plots: dict[str, str] = {}
    if "composite_score" in df.columns:
        fig, ax = plt.subplots(figsize=(max(6, len(df) * 0.5), 4))
        ax.barh(df.index, df["composite_score"], color="#3498db")
        ax.set_xlabel("Composite Score")
        ax.set_title("Segmentation Quality per ROI")
        fig.tight_layout()
        plots["seg_composite_bar"] = fig_to_base64(fig)

    return {"df": df, "plots": plots}


def compute_integration_section(
    adata: AnnData,
    embedding_keys: dict[str, str],
    batch_key: str,
    celltype_key: str | None = None,
    comparison_df: pd.DataFrame | None = None,
) -> dict | None:
    """Compute integration comparison and generate plots.

    Parameters
    ----------
    adata
        AnnData with embeddings.
    embedding_keys
        Mapping of method name to ``obsm`` key.
    batch_key
        Column in ``adata.obs`` with batch labels.
    celltype_key
        Column in ``adata.obs`` with cell type labels. ``None`` to skip bio
        metrics.
    comparison_df
        Pre-computed benchmark DataFrame. When provided, skips calling
        ``compare_integrations()`` and uses this directly for plotting.

    Returns
    -------
    dict or None
        Dict with ``"comparison_df"`` and ``"plots"`` (dict of name to HTML
        or base64 strings). Returns ``None`` on failure.
    """
    if not embedding_keys:
        logger.info("No embeddings provided for integration comparison")
        return None

    try:
        from sc_tools.pl.benchmarking import (
            plot_integration_comparison_table,
            plot_integration_radar,
            plot_integration_ranking_bar,
        )
    except ImportError:
        logger.warning("Integration benchmarking modules not available")
        return None

    if comparison_df is None:
        try:
            from sc_tools.bm.integration import compare_integrations

            comparison_df = compare_integrations(
                adata,
                embedding_keys,
                batch_key=batch_key,
                celltype_key=celltype_key,
                include_unintegrated=True,
                use_scib="sklearn",
            )
        except ImportError:
            logger.warning("sc_tools.bm.integration not available")
            return None
        except Exception:
            logger.warning("Integration comparison failed", exc_info=True)
            return None

    if comparison_df.empty:
        return None

    plots: dict[str, str] = {}

    # Plotly-based interactive plots
    try:
        fig_table = plot_integration_comparison_table(comparison_df)
        plots["integration_table"] = plotly_to_html(fig_table)
    except Exception:
        logger.debug("Integration table plot failed", exc_info=True)

    try:
        fig_radar = plot_integration_radar(comparison_df)
        plots["integration_radar"] = plotly_to_html(fig_radar)
    except Exception:
        logger.debug("Integration radar plot failed", exc_info=True)

    try:
        fig_rank = plot_integration_ranking_bar(comparison_df)
        plots["integration_ranking"] = plotly_to_html(fig_rank)
    except Exception:
        logger.debug("Integration ranking plot failed", exc_info=True)

    # Batch vs bio scatter (if both scores present)
    if "batch_score" in comparison_df.columns and "bio_score" in comparison_df.columns:
        try:
            from sc_tools.pl.benchmarking import plot_batch_vs_bio

            fig_bvb = plot_batch_vs_bio(comparison_df)
            plots["batch_vs_bio"] = plotly_to_html(fig_bvb)
        except Exception:
            logger.debug("Batch vs bio plot failed", exc_info=True)

    return {"comparison_df": comparison_df, "plots": plots}
