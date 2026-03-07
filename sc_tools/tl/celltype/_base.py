"""
Base protocol, registry, and storage helpers for automated cell typing.

Backends must satisfy the CelltypeBackend Protocol.
Results are written to:
- adata.obs[f'celltype_auto_{method}']  (Categorical labels)
- adata.obs[f'celltype_auto_{method}_score']  (float64 confidence)
- adata.obsm[f'celltype_proba_{method}']  (DataFrame, when store_proba=True)
- adata.uns[f'celltype_auto_{method}']  (dict with method, date, metadata)
"""

from __future__ import annotations

from datetime import datetime
from typing import TYPE_CHECKING, Protocol, runtime_checkable

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import anndata as ad


@runtime_checkable
class CelltypeBackend(Protocol):
    """Protocol every cell-typing backend must satisfy."""

    @staticmethod
    def run(
        adata: ad.AnnData,
        *,
        cluster_key: str,
        store_proba: bool,
        **kwargs,
    ) -> tuple[pd.Series, pd.Series, pd.DataFrame | None, dict]:
        """
        Run cell-type annotation.

        Parameters
        ----------
        adata
            Input AnnData (read-only; do not modify in-place).
        cluster_key
            Key in adata.obs holding cluster assignments.
        store_proba
            Whether the caller wants per-cell probability scores.
        **kwargs
            Backend-specific arguments.

        Returns
        -------
        labels
            Cell-level string labels, indexed by obs_names.
        scores
            Cell-level float64 confidence scores, indexed by obs_names.
        proba
            Per-cell per-celltype DataFrame (or None).
        metadata
            Arbitrary dict merged into uns entry.
        """
        ...


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

_BACKENDS: dict[str, type] = {}


def register_celltype_backend(name: str, cls: type) -> None:
    """Register a backend class under a short name."""
    _BACKENDS[name] = cls


def get_backend(name: str) -> type:
    """Return registered backend class or raise ValueError."""
    if name not in _BACKENDS:
        known = sorted(_BACKENDS)
        raise ValueError(
            f"Unknown cell-typing method '{name}'. Available: {known}. "
            "Register custom backends with register_celltype_backend()."
        )
    return _BACKENDS[name]


def list_celltype_methods() -> list[str]:
    """Return alphabetically sorted list of registered backend names."""
    return sorted(_BACKENDS)


# ---------------------------------------------------------------------------
# Storage helper
# ---------------------------------------------------------------------------


def _store_results(
    adata: ad.AnnData,
    method: str,
    labels: pd.Series,
    scores: pd.Series,
    proba: pd.DataFrame | None,
    metadata: dict,
    *,
    store_proba: bool,
) -> None:
    """Write cell-typing results into adata (in-place)."""
    obs_key = f"celltype_auto_{method}"
    score_key = f"celltype_auto_{method}_score"
    obsm_key = f"celltype_proba_{method}"
    uns_key = f"celltype_auto_{method}"

    # Labels as Categorical
    adata.obs[obs_key] = pd.Categorical(labels.reindex(adata.obs_names))

    # Scores as float64
    adata.obs[score_key] = scores.reindex(adata.obs_names).astype(np.float64).values

    # Optional probability matrix
    if store_proba and proba is not None:
        adata.obsm[obsm_key] = proba.reindex(adata.obs_names).values

    # Metadata entry
    uns_entry = {"method": method, "date": datetime.now().isoformat()}
    uns_entry.update(metadata)
    adata.uns[uns_key] = uns_entry
