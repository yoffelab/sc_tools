"""Multi-omic assembly module.

Public API for building MultiOmicAtlas from per-modality AnnData objects.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

__all__ = [
    "MultiOmicAtlas",
    "build_mudata",
    "join_subject_metadata",
    "celltype_proportions",
    "aggregate_by_level",
    "get_embedding_backend",
    "list_embedding_methods",
]


def __getattr__(name: str):
    """Lazy imports for heavy dependencies."""
    if name == "MultiOmicAtlas":
        from sc_tools.assembly._atlas import MultiOmicAtlas

        return MultiOmicAtlas
    if name == "build_mudata":
        from sc_tools.assembly._build import build_mudata

        return build_mudata
    if name == "join_subject_metadata":
        from sc_tools.assembly._metadata import join_subject_metadata

        return join_subject_metadata
    if name == "celltype_proportions":
        from sc_tools.assembly._query import celltype_proportions

        return celltype_proportions
    if name == "aggregate_by_level":
        from sc_tools.assembly._query import aggregate_by_level

        return aggregate_by_level
    if name == "get_embedding_backend":
        from sc_tools.assembly.embed._base import get_embedding_backend

        return get_embedding_backend
    if name == "list_embedding_methods":
        from sc_tools.assembly.embed._base import list_embedding_methods

        return list_embedding_methods
    raise AttributeError(f"module 'sc_tools.assembly' has no attribute {name!r}")
