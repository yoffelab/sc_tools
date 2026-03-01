"""Normalization, transformation, and gene filtering for preprocessing.

Provides:
- normalize_total: Library-size normalization (scanpy / rapids-singlecell).
- log_transform: Log1p transformation.
- scale: Zero-center and scale features.
- arcsinh_transform: Arcsinh transform for mass cytometry (IMC) protein data.
- filter_genes_by_pattern: Remove genes matching regex patterns (MT, ribosomal, hemoglobin).
- backup_raw: Save a copy of adata to adata.raw before any transformation.
"""

from __future__ import annotations

import logging
import re
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from ._gpu import get_backend

logger = logging.getLogger(__name__)

__all__ = [
    "normalize_total",
    "log_transform",
    "scale",
    "arcsinh_transform",
    "filter_genes_by_pattern",
    "backup_raw",
]

# Default gene patterns to exclude: mitochondrial, ribosomal, hemoglobin
DEFAULT_FILTER_PATTERNS = [
    r"^MT-",  # mitochondrial
    r"^RP[SL]",  # ribosomal
    r"^HB[^(P)]",  # hemoglobin (but not HBEGF, HBP1, etc.)
]


def backup_raw(adata: AnnData) -> None:
    """Save a copy of the current adata to adata.raw (no-op if already set).

    Parameters
    ----------
    adata
        Annotated data matrix. Modified in place.
    """
    if adata.raw is not None:
        logger.info("adata.raw already exists; skipping backup")
        return
    adata.raw = adata.copy()
    logger.info("Backed up adata to adata.raw (%d cells x %d genes)", adata.n_obs, adata.n_vars)


def normalize_total(
    adata: AnnData,
    target_sum: float | None = 1e4,
    inplace: bool = True,
    **kwargs: Any,
) -> AnnData | None:
    """Library-size normalize counts per cell.

    Wraps ``scanpy.pp.normalize_total`` (or ``rapids_singlecell.pp.normalize_total``
    on GPU).

    Parameters
    ----------
    adata
        Annotated data with raw counts in ``X``.
    target_sum
        Target total counts per cell after normalization.
    inplace
        If True, modify adata in place.
    **kwargs
        Passed to the backend ``normalize_total``.

    Returns
    -------
    AnnData or None
        Modified adata if ``inplace=False``, else None.
    """
    backend, name = get_backend()
    logger.info("normalize_total (target_sum=%s, backend=%s)", target_sum, name)
    return backend.pp.normalize_total(adata, target_sum=target_sum, inplace=inplace, **kwargs)


def log_transform(
    adata: AnnData,
    base: float | None = None,
    inplace: bool = True,
    **kwargs: Any,
) -> AnnData | None:
    """Apply log1p transformation.

    Wraps ``scanpy.pp.log1p`` (or ``rapids_singlecell.pp.log1p`` on GPU).

    Parameters
    ----------
    adata
        Annotated data (typically after ``normalize_total``).
    base
        Logarithm base. None for natural log (default).
    inplace
        If True, modify adata in place.
    **kwargs
        Passed to the backend ``log1p``.
    """
    backend, name = get_backend()
    logger.info("log1p transform (backend=%s)", name)
    return backend.pp.log1p(adata, base=base, **kwargs)


def scale(
    adata: AnnData,
    max_value: float | None = 10,
    zero_center: bool = True,
    inplace: bool = True,
    **kwargs: Any,
) -> AnnData | None:
    """Zero-center and scale features to unit variance.

    Wraps ``scanpy.pp.scale`` (or ``rapids_singlecell.pp.scale`` on GPU).

    Parameters
    ----------
    adata
        Annotated data.
    max_value
        Clip values to this maximum after scaling (default 10).
    zero_center
        If True, center each gene to zero mean.
    inplace
        If True, modify adata in place.
    **kwargs
        Passed to the backend ``scale``.
    """
    backend, name = get_backend()
    logger.info("scale (max_value=%s, backend=%s)", max_value, name)
    return backend.pp.scale(adata, max_value=max_value, zero_center=zero_center, **kwargs)


def arcsinh_transform(
    adata: AnnData,
    cofactor: float = 5,
    inplace: bool = True,
) -> AnnData | None:
    """Arcsinh transform for mass cytometry (IMC) protein data.

    Applies ``arcsinh(X / cofactor)`` element-wise. This is the standard
    normalization for CyTOF / IMC data (NOT log1p).

    Parameters
    ----------
    adata
        Annotated data with raw protein intensities in ``X``.
    cofactor
        Scaling factor before arcsinh. Standard value is 5 for CyTOF/IMC.
    inplace
        If True, modify adata.X in place.

    Returns
    -------
    AnnData or None
        Modified adata if ``inplace=False``, else None.
    """
    logger.info("arcsinh transform (cofactor=%s)", cofactor)

    if not inplace:
        adata = adata.copy()

    if sparse.issparse(adata.X):
        adata.X = np.arcsinh(adata.X.toarray() / cofactor)
    else:
        adata.X = np.arcsinh(adata.X / cofactor)

    if not inplace:
        return adata
    return None


def filter_genes_by_pattern(
    adata: AnnData,
    patterns: list[str] | None = None,
    exclude: bool = True,
    case_sensitive: bool = False,
) -> None:
    """Remove (or keep) genes matching regex patterns in place.

    Parameters
    ----------
    adata
        Annotated data. Modified in place.
    patterns
        List of regex patterns. Defaults to ``["^MT-", "^RP[SL]", "^HB[^(P)]"]``
        (mitochondrial, ribosomal, hemoglobin).
    exclude
        If True (default), remove matching genes. If False, keep only matching genes.
    case_sensitive
        If False (default), patterns are case-insensitive.
    """
    if patterns is None:
        patterns = DEFAULT_FILTER_PATTERNS

    flags = 0 if case_sensitive else re.IGNORECASE
    var_names = pd.Series(adata.var_names)

    mask = pd.Series(False, index=var_names.index)
    for pattern in patterns:
        mask |= var_names.str.contains(pattern, regex=True, flags=flags)

    n_matching = mask.sum()
    if exclude:
        keep = ~mask
        logger.info(
            "Removing %d/%d genes matching patterns %s",
            n_matching,
            adata.n_vars,
            patterns,
        )
    else:
        keep = mask
        logger.info(
            "Keeping %d/%d genes matching patterns %s",
            n_matching,
            adata.n_vars,
            patterns,
        )

    adata._inplace_subset_var(keep.values)
