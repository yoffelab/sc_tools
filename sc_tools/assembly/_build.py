"""MuData construction from per-modality AnnData objects.

Provides:
- build_mudata: Create a MuData object from a dict of AnnData objects,
  with patient and sample metadata stored in uns.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData

__all__ = ["build_mudata"]


def build_mudata(
    modalities: dict[str, AnnData],
    *,
    subject_key: str = "subject_id",
):
    """Build a MuData object from per-modality AnnData objects.

    Parameters
    ----------
    modalities
        Dict mapping modality name to AnnData.
    subject_key
        Column name for subject identifiers.

    Returns
    -------
    mudata.MuData
        Multi-modal data container with patient/sample metadata in uns.
    """
    import mudata  # lazy import per CLI-08

    from sc_tools.assembly._metadata import join_subject_metadata

    # Validate and build metadata
    patient_metadata, sample_metadata = join_subject_metadata(
        modalities, subject_key=subject_key
    )

    # Create MuData -- keys are modality names per D-04
    mdata = mudata.MuData(modalities)

    # Sync global obs (Pitfall 1)
    mdata.update()

    # Store metadata in uns
    mdata.uns["patient_metadata"] = patient_metadata
    mdata.uns["sample_metadata"] = sample_metadata

    return mdata
