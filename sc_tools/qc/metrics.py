"""
QC metrics: wrappers around scanpy for spots/cells and genes.

Provides:
- calculate_qc_metrics: total_counts, n_genes_by_counts, pct_counts_mt (and optionally hb).
- filter_cells, filter_genes: count-based filtering.
- highly_variable_genes: HVG selection (wraps scanpy).

Two usage points: pre-normalization (raw counts) and post-normalization (optional second pass).
"""

from __future__ import annotations

import re
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData

__all__ = [
    "calculate_qc_metrics",
    "filter_cells",
    "filter_genes",
    "highly_variable_genes",
]


def calculate_qc_metrics(
    adata: AnnData,
    *,
    mt_pattern: str | re.Pattern | None = "^(MT-|mt-|Mt)",
    hb_pattern: str | re.Pattern | None = "^(HB|Hb|HBEGF)",
    qc_vars: list[str] | None = None,
    inplace: bool = True,
    percent_top: list[int] | None = (50, 100, 200, 500),
    log1p: bool = False,
    modality: str = "visium",
    **kwargs: Any,
) -> pd.DataFrame | None:
    """
    Compute QC metrics for spots/cells and genes (wrap scanpy).

    Optionally marks mitochondrial and hemoglobin genes and adds
    pct_counts_mt (and pct_counts_hb) to adata.obs.

    Parameters
    ----------
    adata : AnnData
        Annotated data (raw counts recommended).
    mt_pattern : str or compiled regex or None
        Pattern to mark mitochondrial genes in adata.var (default MT- / mt- / Mt).
        If None, mt genes are not marked and pct_counts_mt is not computed.
        Automatically set to None for protein-based modalities (e.g. ``"imc"``).
    hb_pattern : str or compiled regex or None
        Pattern to mark hemoglobin genes for pct_counts_hb. If None, not computed.
        Automatically set to None for protein-based modalities.
    qc_vars : list of str or None
        If None, built from mt_pattern and hb_pattern: ['mt'] and optionally ['mt','hb'].
        Passed to scanpy as qc_vars (column names in adata.var).
    inplace : bool
        If True, add metrics to adata.obs and adata.var (default True).
    percent_top : list of int or None
        Passed to scanpy (default (50, 100, 200, 500)).
        Automatically capped to ``n_vars`` to avoid IndexError on small panels.
    log1p : bool
        Passed to scanpy (default False).
    modality : str
        Data modality. Protein-based modalities (``"imc"``) skip MT/HB patterns.
    **kwargs
        Passed to scanpy.pp.calculate_qc_metrics.

    Returns
    -------
    DataFrame or None
        If inplace is False, returns (obs_df, var_df) concatenation as per scanpy;
        otherwise None.
    """
    import scanpy as sc

    # Protein-based modalities have no mitochondrial or hemoglobin genes
    _protein_modalities = {"imc"}
    if modality in _protein_modalities:
        mt_pattern = None
        hb_pattern = None

    var_names = pd.Series(adata.var_names)
    if mt_pattern is not None:
        if isinstance(mt_pattern, str):
            mt_pattern = re.compile(mt_pattern, re.IGNORECASE)
        adata.var["mt"] = var_names.str.match(mt_pattern).values
    if hb_pattern is not None:
        if isinstance(hb_pattern, str):
            hb_pattern = re.compile(hb_pattern, re.IGNORECASE)
        adata.var["hb"] = var_names.str.match(hb_pattern).values

    if qc_vars is None:
        qc_vars = []
        if mt_pattern is not None and "mt" in adata.var.columns:
            qc_vars.append("mt")
        if hb_pattern is not None and "hb" in adata.var.columns:
            qc_vars.append("hb")

    # Cap percent_top to n_vars to avoid IndexError on small panels (e.g. 52-protein IMC)
    if percent_top is not None:
        n_vars = adata.n_vars
        percent_top = tuple(p for p in percent_top if p <= n_vars) or None

    return sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=qc_vars,
        inplace=inplace,
        percent_top=percent_top,
        log1p=log1p,
        **kwargs,
    )


def filter_cells(
    adata: AnnData,
    min_counts: int | None = None,
    min_genes: int | None = None,
    max_counts: int | None = None,
    max_genes: int | None = None,
    inplace: bool = True,
    **kwargs: Any,
) -> AnnData | None:
    """
    Filter out spots/cells by counts and number of genes (wrap scanpy).

    Parameters
    ----------
    adata : AnnData
        Annotated data (should have total_counts and n_genes_by_counts in obs
        from calculate_qc_metrics).
    min_counts, min_genes, max_counts, max_genes : int or None
        Thresholds; None means do not apply. Passed to scanpy.
    inplace : bool
        If True, filter in place (default True).
    **kwargs
        Passed to scanpy.pp.filter_cells.

    Returns
    -------
    AnnData or None
        Filtered object if inplace=False, else None.
    """
    import scanpy as sc

    return sc.pp.filter_cells(
        adata,
        min_counts=min_counts,
        min_genes=min_genes,
        max_counts=max_counts,
        max_genes=max_genes,
        inplace=inplace,
        **kwargs,
    )


def filter_genes(
    adata: AnnData,
    min_counts: int | None = None,
    min_cells: int | None = None,
    max_counts: int | None = None,
    max_cells: int | None = None,
    inplace: bool = True,
    **kwargs: Any,
) -> AnnData | None:
    """
    Filter out genes by counts and number of cells (wrap scanpy).

    Parameters
    ----------
    adata : AnnData
        Annotated data.
    min_counts, min_cells, max_counts, max_cells : int or None
        Thresholds; None means do not apply. Passed to scanpy.
    inplace : bool
        If True, filter in place (default True).
    **kwargs
        Passed to scanpy.pp.filter_genes.

    Returns
    -------
    AnnData or None
        Filtered object if inplace=False, else None.
    """
    import scanpy as sc

    return sc.pp.filter_genes(
        adata,
        min_counts=min_counts,
        min_cells=min_cells,
        max_counts=max_counts,
        max_cells=max_cells,
        inplace=inplace,
        **kwargs,
    )


def highly_variable_genes(
    adata: AnnData,
    flavor: str = "seurat",
    n_top_genes: int | None = None,
    min_mean: float = 0.0125,
    max_mean: float = 3,
    min_disp: float = 0.5,
    max_disp: float = np.inf,
    batch_key: str | None = None,
    subset: bool = False,
    inplace: bool = True,
    **kwargs: Any,
) -> pd.DataFrame | None:
    """
    Mark or subset highly variable genes (wrap scanpy).

    Parameters
    ----------
    adata : AnnData
        Annotated data (normalized, e.g. log1p, recommended for seurat_v3).
    flavor : str
        'seurat', 'seurat_v3', or 'cell_ranger' (default 'seurat').
    n_top_genes : int or None
        If set, use this many top genes (common for seurat_v3).
    min_mean, max_mean, min_disp, max_disp : float
        Passed to scanpy (flavor-dependent).
    batch_key : str or None
        If set, compute HVGs per batch (e.g. 'sample').
    subset : bool
        If True, subset adata to HVGs (default False).
    inplace : bool
        If True, add 'highly_variable' etc. to adata.var (default True).
    **kwargs
        Passed to scanpy.pp.highly_variable_genes.

    Returns
    -------
    DataFrame or None
        If inplace=False, returns the var DataFrame with HVG columns; else None.
    """
    import scanpy as sc

    return sc.pp.highly_variable_genes(
        adata,
        flavor=flavor,
        n_top_genes=n_top_genes,
        min_mean=min_mean,
        max_mean=max_mean,
        min_disp=min_disp,
        max_disp=max_disp,
        batch_key=batch_key,
        subset=subset,
        inplace=inplace,
        **kwargs,
    )
