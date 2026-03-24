"""
Subject-level metadata validation for multi-sample projects.

Provides:
- validate_subject_metadata: Check subject_id presence, distinctness, confounding.
- check_confounding: Detect perfect batch-condition confounding via crosstab.

These functions return warning strings (never raise) per D-08, D-09, D-10.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from anndata import AnnData

__all__ = ["validate_subject_metadata", "check_confounding"]

logger = logging.getLogger(__name__)


def validate_subject_metadata(
    adata: AnnData,
    *,
    multi_sample: bool | None = None,
    subject_key: str = "subject_id",
    library_key: str = "library_id",
    condition_key: str | None = None,
    batch_key: str | None = None,
) -> list[str]:
    """Validate subject-level metadata in an AnnData object.

    Parameters
    ----------
    adata
        Annotated data matrix.
    multi_sample
        Whether this is a multi-sample project. If None, auto-detected from
        the number of unique library_id values.
    subject_key
        Column name for subject identifiers.
    library_key
        Column name for library identifiers.
    condition_key
        Column name for experimental condition (for confounding check).
    batch_key
        Column name for batch (for confounding check).

    Returns
    -------
    list[str]
        Warning messages. Empty list means all validations passed.
    """
    warnings: list[str] = []

    # Auto-detect multi-sample
    if multi_sample is None:
        multi_sample = (
            library_key in adata.obs.columns
            and adata.obs[library_key].nunique() > 1
        )

    if multi_sample:
        if subject_key not in adata.obs.columns:
            warnings.append(
                f"Multi-sample project but '{subject_key}' not found in obs. "
                "subject_id is required for pseudobulk DE and prevents pseudoreplication."
            )
        elif library_key in adata.obs.columns:
            # Check distinctness: subject_key should not map 1:1 to library_key
            mapping = (
                adata.obs[[subject_key, library_key]]
                .drop_duplicates()
            )
            if mapping.shape[0] == adata.obs[library_key].nunique():
                # Every library maps to exactly one subject and vice versa
                subject_nunique = adata.obs[subject_key].nunique()
                library_nunique = adata.obs[library_key].nunique()
                if subject_nunique == library_nunique:
                    warnings.append(
                        f"'{subject_key}' maps 1:1 to '{library_key}' -- "
                        "these may be the same. subject_id should represent "
                        "biological replicates."
                    )
    else:
        if subject_key not in adata.obs.columns:
            warnings.append(
                f"'{subject_key}' not found in obs "
                "(single-sample: informational only)."
            )

    # Batch-condition confounding check
    if (
        condition_key
        and batch_key
        and condition_key in adata.obs.columns
        and batch_key in adata.obs.columns
    ):
        if check_confounding(adata.obs, batch_key, condition_key):
            warnings.append(
                "Perfect batch-condition confounding detected: every batch "
                "maps to exactly one condition level. DE results may be unreliable."
            )

    return warnings


def check_confounding(
    obs: pd.DataFrame,
    batch_key: str,
    condition_key: str,
) -> bool:
    """Check if batch perfectly confounds condition.

    Returns True if every batch maps to exactly one condition level,
    meaning batch and condition effects cannot be separated.

    Parameters
    ----------
    obs
        DataFrame with batch and condition columns.
    batch_key
        Column name for batch.
    condition_key
        Column name for condition.

    Returns
    -------
    bool
        True if perfect confounding detected.
    """
    crosstab = pd.crosstab(obs[batch_key], obs[condition_key])
    # Perfect confounding: each batch has cells from exactly one condition
    return bool((crosstab > 0).sum(axis=1).eq(1).all())
