"""
Marker validation for post-celltyping QC.

Checks whether assigned cell types express their canonical marker genes
at expected levels. Types where ALL canonical markers fall below a
threshold are flagged as potentially mis-annotated.

Flagging is informational only -- this module never raises on low expression.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from anndata import AnnData

__all__ = ["compute_marker_validation", "render_marker_dotplot"]

logger = logging.getLogger(__name__)


def compute_marker_validation(
    adata: AnnData,
    celltype_key: str,
    marker_genes: dict[str, list[str]],
    *,
    threshold: float = 0.1,
    top_n: int = 5,
) -> tuple[pd.DataFrame, dict]:
    """Validate cell type assignments against canonical marker expression.

    For each cell type in *marker_genes*, computes mean expression of each
    marker gene across cells of that type. A cell type is flagged when
    **all** its canonical markers have mean expression below *threshold*.

    Parameters
    ----------
    adata
        Annotated data matrix.
    celltype_key
        Column in ``adata.obs`` with cell type labels.
    marker_genes
        Dict mapping cell type name to list of marker gene symbols.
    threshold
        Expression threshold below which a marker is considered low.
        A cell type is flagged only if ALL its markers are below this.
    top_n
        Maximum number of markers per type to include (unused in compute,
        reserved for dotplot rendering).

    Returns
    -------
    tuple[pd.DataFrame, dict]
        - DataFrame with columns: ``celltype``, ``marker_gene``,
          ``mean_expr``, ``flagged``.
        - Summary dict with keys: ``n_types_tested``, ``n_flagged``,
          ``total_cells``.
    """
    if not marker_genes:
        return pd.DataFrame(columns=["celltype", "marker_gene", "mean_expr", "flagged"]), {
            "n_types_tested": 0,
            "n_flagged": 0,
            "total_cells": int(adata.n_obs),
        }

    rows: list[dict] = []
    flagged_types: set[str] = set()
    types_tested: set[str] = set()

    for ct, genes in marker_genes.items():
        # Only test types that exist in the data
        if ct not in adata.obs[celltype_key].cat.categories if hasattr(adata.obs[celltype_key], "cat") else ct not in adata.obs[celltype_key].unique():
            continue

        types_tested.add(ct)
        mask = adata.obs[celltype_key] == ct
        subset = adata[mask]

        all_below = True
        for gene in genes:
            if gene in adata.var_names:
                gene_idx = list(adata.var_names).index(gene)
                expr_values = subset.X[:, gene_idx]
                # Handle sparse matrices
                if hasattr(expr_values, "toarray"):
                    expr_values = expr_values.toarray().flatten()
                elif hasattr(expr_values, "A"):
                    expr_values = expr_values.A.flatten()
                else:
                    expr_values = np.asarray(expr_values).flatten()
                mean_expr = float(np.nanmean(expr_values))
                if mean_expr >= threshold:
                    all_below = False
            else:
                mean_expr = float("nan")
                # NaN genes don't count toward the "all below" decision

            rows.append({
                "celltype": ct,
                "marker_gene": gene,
                "mean_expr": mean_expr,
                "flagged": False,  # placeholder, set after loop
            })

        if all_below:
            flagged_types.add(ct)

    df = pd.DataFrame(rows)

    # Set flagged column based on whether the celltype is flagged
    if len(df) > 0:
        df["flagged"] = df["celltype"].isin(flagged_types)

    summary = {
        "n_types_tested": len(types_tested),
        "n_flagged": len(flagged_types),
        "total_cells": int(adata.n_obs),
    }

    return df, summary


def render_marker_dotplot(
    adata: AnnData,
    celltype_key: str,
    marker_genes: dict[str, list[str]],
    *,
    top_n: int = 5,
) -> str:
    """Render a dotplot of marker gene expression per cell type.

    Uses ``scanpy.pl.dotplot`` to render top *top_n* markers per cell type,
    returning the result as a base64-encoded PNG string.

    Parameters
    ----------
    adata
        Annotated data matrix.
    celltype_key
        Column in ``adata.obs`` with cell type labels.
    marker_genes
        Dict mapping cell type name to list of marker gene symbols.
    top_n
        Maximum number of markers per cell type to display.

    Returns
    -------
    str
        Base64-encoded PNG string, or empty string if scanpy is unavailable
        or the plot fails.
    """
    try:
        import scanpy as sc

        from .report_utils import fig_to_base64
    except ImportError:
        logger.debug("scanpy not available for marker dotplot")
        return ""

    # Trim to top_n markers per type and filter to genes in var_names
    trimmed: dict[str, list[str]] = {}
    for ct, genes in marker_genes.items():
        valid = [g for g in genes[:top_n] if g in adata.var_names]
        if valid:
            trimmed[ct] = valid

    if not trimmed:
        return ""

    try:
        fig_dot, _ = sc.pl.dotplot(
            adata,
            var_names=trimmed,
            groupby=celltype_key,
            return_fig=True,
            show=False,
        )
        return fig_to_base64(fig_dot)
    except Exception:
        logger.debug("Marker dotplot rendering failed", exc_info=True)
        return ""
