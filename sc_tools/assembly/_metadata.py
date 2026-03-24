"""Hierarchical metadata join for multi-omic assembly.

Provides:
- join_subject_metadata: Outer join of subject_ids across modalities,
  producing patient_metadata and sample_metadata DataFrames.
"""

from __future__ import annotations

import logging
import warnings
from itertools import combinations
from typing import TYPE_CHECKING

import pandas as pd

from sc_tools.errors import SCToolsDataError

if TYPE_CHECKING:
    from anndata import AnnData

__all__ = ["join_subject_metadata"]

logger = logging.getLogger(__name__)


def join_subject_metadata(
    modalities: dict[str, AnnData],
    *,
    subject_key: str = "subject_id",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build patient and sample metadata via outer join across modalities.

    Parameters
    ----------
    modalities
        Dict mapping modality name to AnnData (e.g. ``{"rna": adata_rna, ...}``).
    subject_key
        Column name for subject identifiers in ``obs``.

    Returns
    -------
    patient_metadata
        DataFrame with one row per unique subject, boolean ``has_{mod}`` columns.
    sample_metadata
        DataFrame with one row per unique (sample_id, modality) pair.

    Raises
    ------
    SCToolsDataError
        If any modality is missing the ``subject_key`` column.
    """
    # Validate subject_key presence
    for mod_name, adata in modalities.items():
        if subject_key not in adata.obs.columns:
            raise SCToolsDataError(
                f"Modality '{mod_name}' is missing '{subject_key}' in obs. "
                f"All modalities must have a '{subject_key}' column for assembly.",
                suggestion=f"Add '{subject_key}' to adata.obs before assembly.",
            )

    # Run validate_subject_metadata on each modality if available
    try:
        from sc_tools.qc.metadata import validate_subject_metadata

        for mod_name, adata in modalities.items():
            warns = validate_subject_metadata(adata, subject_key=subject_key)
            for w in warns:
                logger.warning("Modality '%s': %s", mod_name, w)
    except ImportError:
        pass

    # Collect unique subject_ids per modality
    mod_subjects: dict[str, set[str]] = {}
    for mod_name, adata in modalities.items():
        mod_subjects[mod_name] = set(adata.obs[subject_key].unique())

    # Check for zero overlap between any pair (Pitfall 4)
    all_subjects = set()
    for subjects in mod_subjects.values():
        all_subjects.update(subjects)

    for mod_a, mod_b in combinations(mod_subjects.keys(), 2):
        overlap = mod_subjects[mod_a] & mod_subjects[mod_b]
        if len(overlap) == 0:
            warnings.warn(
                f"No overlap in {subject_key} between modalities "
                f"'{mod_a}' and '{mod_b}'. Outer join will include all patients "
                "but cross-modal queries will return empty for non-shared patients.",
                UserWarning,
                stacklevel=2,
            )

    # Build patient_metadata: outer join (union of all subjects)
    patient_rows = []
    for subj in sorted(all_subjects):
        row: dict[str, object] = {subject_key: subj}
        for mod_name in modalities:
            row[f"has_{mod_name}"] = subj in mod_subjects[mod_name]
        patient_rows.append(row)
    patient_metadata = pd.DataFrame(patient_rows)

    # Build sample_metadata: concatenation of per-modality sample info
    sample_rows = []
    for mod_name, adata in modalities.items():
        if "sample_id" in adata.obs.columns:
            sample_info = (
                adata.obs[[subject_key, "sample_id"]]
                .drop_duplicates()
                .assign(modality=mod_name)
            )
            sample_rows.append(sample_info)
        else:
            # If no sample_id, use subject_id as sample_id
            sample_info = (
                adata.obs[[subject_key]]
                .drop_duplicates()
                .assign(sample_id=lambda df: df[subject_key], modality=mod_name)
            )
            sample_rows.append(sample_info)

    sample_metadata = pd.concat(sample_rows, ignore_index=True) if sample_rows else pd.DataFrame()

    return patient_metadata, sample_metadata
