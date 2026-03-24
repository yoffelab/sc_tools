"""Embedding backend Protocol, registry, and dispatch.

Mirrors the celltype backend pattern from sc_tools.tl.celltype._base.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Protocol, runtime_checkable

import numpy as np

if TYPE_CHECKING:
    pass

__all__ = [
    "EmbeddingBackend",
    "register_embedding_backend",
    "get_embedding_backend",
    "list_embedding_methods",
]


@runtime_checkable
class EmbeddingBackend(Protocol):
    """Protocol every joint embedding backend must satisfy."""

    @staticmethod
    def run(mdata, *, n_factors: int = 15, **kwargs) -> tuple[np.ndarray, dict]:
        """Run joint embedding.

        Parameters
        ----------
        mdata
            MuData object with multi-modal data.
        n_factors
            Number of latent factors/dimensions.
        **kwargs
            Backend-specific arguments.

        Returns
        -------
        embedding
            Array of shape (n_cells, n_factors).
        metadata
            Dict with method name, parameters, etc.
        """
        ...


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

_BACKENDS: dict[str, type] = {}


def register_embedding_backend(name: str, cls: type) -> None:
    """Register a backend class under a short name."""
    _BACKENDS[name] = cls


def get_embedding_backend(name: str) -> type:
    """Return registered backend class or raise ValueError."""
    if name not in _BACKENDS:
        known = sorted(_BACKENDS)
        raise ValueError(
            f"Unknown embedding method '{name}'. Available: {known}. "
            "Register custom backends with register_embedding_backend()."
        )
    return _BACKENDS[name]


def list_embedding_methods() -> list[str]:
    """Return alphabetically sorted list of registered backend names."""
    return sorted(_BACKENDS)
