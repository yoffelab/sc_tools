"""Joint embedding backends for multi-omic data.

Exports the embedding backend Protocol, registry functions, and
triggers registration of all built-in backends.
"""

from __future__ import annotations

from sc_tools.assembly.embed._base import (
    EmbeddingBackend,
    get_embedding_backend,
    list_embedding_methods,
    register_embedding_backend,
)

# Import backends to trigger registration
from sc_tools.assembly.embed._mofa import MofaBackend  # noqa: F401
from sc_tools.assembly.embed._multivi import MultiviBackend  # noqa: F401
from sc_tools.assembly.embed._totalvi import TotalviBackend  # noqa: F401

__all__ = [
    "EmbeddingBackend",
    "MofaBackend",
    "MultiviBackend",
    "TotalviBackend",
    "get_embedding_backend",
    "list_embedding_methods",
    "register_embedding_backend",
]
