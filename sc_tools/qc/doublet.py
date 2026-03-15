"""
Doublet detection for single-cell resolution modalities using scVI Solo.

Provides:
- score_doublets_solo: Train Solo on raw counts and annotate obs with doublet scores.

Applicable modalities: xenium, cosmx, visium_hd_cell (single-cell resolution).
NOT applicable: visium (spot-based, ~10-50 cells per spot),
                visium_hd (bin-level, multiple cells per bin), imc.
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
from anndata import AnnData

try:
    from scvi.external import SOLO
    from scvi.model import SCVI
except ImportError:  # pragma: no cover
    SCVI = None  # type: ignore[assignment,misc]
    SOLO = None  # type: ignore[assignment,misc]

__all__ = ["score_doublets_solo"]

logger = logging.getLogger(__name__)

# Modalities for which doublet detection does not apply (not single-cell resolution)
_SPOT_BASED_MODALITIES = frozenset({"visium", "visium_hd", "imc"})

# Modalities that are single-cell resolution and where doublet detection applies
_SC_RESOLUTION_MODALITIES = frozenset({"xenium", "cosmx", "visium_hd_cell"})


def _check_modality(adata: AnnData, modality: str | None) -> None:
    """
    Raise ValueError if modality is not single-cell resolution.

    Modality is determined in priority order:
    1. Explicit ``modality`` parameter (if not None).
    2. ``adata.obs['modality']`` column (if present and all values agree).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    modality : str or None
        Explicit modality override.  If None, inferred from obs.

    Raises
    ------
    ValueError
        If the resolved modality is in the spot-based / non-SC-resolution set.
    """
    resolved: str | None = modality

    if resolved is None and "modality" in adata.obs.columns:
        unique_modalities = adata.obs["modality"].unique().tolist()
        if len(unique_modalities) == 1:
            resolved = str(unique_modalities[0])
        elif len(unique_modalities) > 1:
            # Mixed modalities in the same AnnData -- check if any is spot-based
            for m in unique_modalities:
                if m in _SPOT_BASED_MODALITIES:
                    raise ValueError(
                        f"score_doublets_solo is not applicable to spot-based or "
                        f"non-single-cell-resolution modality '{m}'. "
                        f"Use only for: {sorted(_SC_RESOLUTION_MODALITIES)}."
                    )
            return  # all modalities are SC-resolution -- proceed

    if resolved is not None and resolved in _SPOT_BASED_MODALITIES:
        raise ValueError(
            f"score_doublets_solo is not applicable to spot-based or "
            f"non-single-cell-resolution modality '{resolved}'. "
            f"Use only for: {sorted(_SC_RESOLUTION_MODALITIES)}. "
            f"Visium and Visium HD spots contain multiple cells; "
            f"IMC uses protein-based imaging without transcript-level doublets."
        )


def score_doublets_solo(
    adata: AnnData,
    *,
    batch_key: str = "library_id",
    max_epochs: int = 400,
    use_gpu: bool = True,
    threshold: float = 0.5,
    modality: str | None = None,
    scvi_train_kwargs: dict[str, Any] | None = None,
    solo_train_kwargs: dict[str, Any] | None = None,
) -> AnnData:
    """
    Score cells as doublets using scVI Solo.

    Trains a minimal scVI model on raw counts (or reuses one if
    ``adata.uns['scvi_model_path']`` is set), then runs Solo doublet
    detection.  The function ONLY annotates -- it does NOT filter cells.
    The caller decides whether to remove doublets based on the scores.

    Applicable modalities
    ---------------------
    - ``xenium``: single-cell resolution targeted transcriptomics.
    - ``cosmx``: single-cell resolution spatial transcriptomics.
    - ``visium_hd_cell``: SpaceRanger 4 cell-segmented Visium HD data.

    NOT applicable
    --------------
    - ``visium``: each spot covers approximately 10-50 cells; doublet
      detection is not meaningful for spot mixtures.
    - ``visium_hd``: 8 um bins still contain multiple cells.
    - ``imc``: protein-based imaging; Solo is not trained on protein data.

    Parameters
    ----------
    adata : AnnData
        Annotated data with raw counts in ``adata.X``.  Must NOT be
        normalized -- scVI requires raw integer counts.
    batch_key : str
        Column in ``adata.obs`` used as the batch variable for scVI.
        Defaults to ``'library_id'``.
    max_epochs : int
        Maximum training epochs for both scVI and Solo models.
        Defaults to 400.  Reduce for quick testing.
    use_gpu : bool
        Whether to use GPU acceleration.  Defaults to True.
        Set to False for CPU-only environments.
    threshold : float
        Probability threshold above which a cell is classified as a
        doublet.  Defaults to 0.5.
    modality : str or None
        Explicit modality label.  If None, inferred from
        ``adata.obs['modality']`` when present.  Used for the modality
        guard -- raises ValueError for spot-based modalities.
    scvi_train_kwargs : dict or None
        Extra keyword arguments forwarded to ``SCVI.train()``.
    solo_train_kwargs : dict or None
        Extra keyword arguments forwarded to ``SOLO.train()``.

    Returns
    -------
    AnnData
        The input ``adata`` with the following additions:

        - ``obs['solo_doublet_score']`` (float64): probability of being a
          doublet as returned by Solo.
        - ``obs['solo_is_doublet']`` (bool): True when
          ``solo_doublet_score > threshold``.
        - ``uns['solo_summary']`` (dict): summary statistics with keys
          ``n_doublets`` (int), ``pct_doublets`` (float), and
          ``threshold`` (float).

    Raises
    ------
    ValueError
        If the resolved modality is not single-cell resolution
        (e.g. visium, visium_hd, imc).
    ImportError
        If ``scvi-tools`` is not installed.

    Examples
    --------
    Score doublets for a Xenium dataset:

    >>> adata = score_doublets_solo(adata, batch_key="library_id")
    >>> adata_clean = adata[~adata.obs["solo_is_doublet"]].copy()
    """
    # Modality guard -- runs before any import check so callers get an
    # informative error even when scvi-tools is not installed.
    _check_modality(adata, modality)

    if SCVI is None or SOLO is None:
        raise ImportError(
            "scvi-tools is required for score_doublets_solo. Install with: pip install scvi-tools"
        )

    scvi_kwargs = scvi_train_kwargs or {}
    s_kwargs = solo_train_kwargs or {}

    logger.info(
        "score_doublets_solo: training scVI on %d cells x %d genes "
        "(batch_key=%r, max_epochs=%d, use_gpu=%s)",
        adata.n_obs,
        adata.n_vars,
        batch_key,
        max_epochs,
        use_gpu,
    )

    # Set up scVI model
    batch_key_arg = batch_key if batch_key in adata.obs.columns else None
    SCVI.setup_anndata(adata, batch_key=batch_key_arg)
    scvi_model = SCVI(adata)
    scvi_model.train(max_epochs=max_epochs, use_gpu=use_gpu, **scvi_kwargs)

    # Train Solo on the scVI model
    logger.info("score_doublets_solo: training Solo doublet detector")
    solo_model = SOLO.from_scvi_model(scvi_model)
    solo_model.train(max_epochs=max_epochs, use_gpu=use_gpu, **s_kwargs)

    # Predict doublet probabilities
    predictions = solo_model.predict(soft=True)
    # predictions is a DataFrame with columns ['doublet', 'singlet']
    doublet_scores = np.array(predictions["doublet"], dtype=np.float64)

    is_doublet = doublet_scores > threshold

    # Write results to obs
    adata.obs["solo_doublet_score"] = doublet_scores
    adata.obs["solo_is_doublet"] = is_doublet.astype(bool)

    n_doublets = int(is_doublet.sum())
    pct_doublets = float(100.0 * n_doublets / adata.n_obs)

    adata.uns["solo_summary"] = {
        "n_doublets": n_doublets,
        "pct_doublets": pct_doublets,
        "threshold": float(threshold),
    }

    logger.info(
        "score_doublets_solo: identified %d doublets (%.1f%%) at threshold=%.2f",
        n_doublets,
        pct_doublets,
        threshold,
    )

    return adata
