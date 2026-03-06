"""Shared utilities for h5ad write compatibility in DLBCL pipeline."""

import logging

import pandas as pd

logger = logging.getLogger(__name__)


def clean_adata_for_h5ad(adata):
    """Clean AnnData obs/var dtypes for h5ad write compatibility.

    Seurat-converted objects often have mixed-type columns (numeric + NaN + string)
    that h5py cannot serialize. This forces all columns to either numeric or string.
    Also removes reserved '_index' column from var/raw.var.
    """
    for col in adata.obs.columns:
        s = adata.obs[col]
        numeric = pd.to_numeric(s, errors="coerce")
        if numeric.notna().sum() == s.notna().sum() and s.notna().any():
            adata.obs[col] = numeric
        else:
            adata.obs[col] = s.astype(str).replace(
                {"nan": "", "None": "", "<NA>": ""}
            )
    for col in adata.var.columns:
        s = adata.var[col]
        numeric = pd.to_numeric(s, errors="coerce")
        if numeric.notna().sum() == s.notna().sum() and s.notna().any():
            adata.var[col] = numeric
        else:
            adata.var[col] = s.astype(str).replace(
                {"nan": "", "None": "", "<NA>": ""}
            )
    # Fix _index reserved column
    if "_index" in adata.var.columns:
        adata.var = adata.var.drop(columns=["_index"])
    if adata.raw is not None and "_index" in adata.raw.var.columns:
        raw_var = adata.raw.var.copy()
        raw_var = raw_var.drop(columns=["_index"])
        from anndata import Raw

        adata._raw = Raw(adata, X=adata.raw.X, var=raw_var, varm=adata.raw.varm)

    return adata
