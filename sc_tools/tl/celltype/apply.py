"""
Apply a manually curated cluster-to-celltype mapping to AnnData.

Reads a JSON file or dict of the form::

    {
        "0": {"celltype": "T cell", "celltype_broad": "T", "color": "#e41a1c"},
        "1": {"celltype": "B cell", "celltype_broad": "B"},
        ...
    }

and writes:
- ``adata.obs[celltype_key]``
- ``adata.obs[celltype_broad_key]``
- ``adata.uns[f'{celltype_key}_colors']``  (list aligned with categories)
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import anndata as ad


# A reasonably large color palette for automatic assignment
_DEFAULT_PALETTE = [
    "#e41a1c",
    "#377eb8",
    "#4daf4a",
    "#984ea3",
    "#ff7f00",
    "#a65628",
    "#f781bf",
    "#999999",
    "#66c2a5",
    "#fc8d62",
    "#8da0cb",
    "#e78ac3",
    "#a6d854",
    "#ffd92f",
    "#e5c494",
    "#b3b3b3",
]


def apply_celltype_map(
    adata: ad.AnnData,
    mapping: dict | str | Path,
    cluster_key: str = "leiden",
    celltype_key: str = "celltype",
    celltype_broad_key: str = "celltype_broad",
    unknown_label: str = "Unknown",
    unknown_color: str = "#cccccc",
    copy: bool = False,
) -> ad.AnnData:
    """
    Apply a cluster-to-celltype mapping to AnnData.

    Parameters
    ----------
    adata
        Annotated data matrix.
    mapping
        Dict or path to JSON file. Keys are cluster IDs (as strings); values
        are dicts with at minimum a ``'celltype'`` field, optionally
        ``'celltype_broad'`` and ``'color'``.
    cluster_key
        Column in ``adata.obs`` with cluster assignments.
    celltype_key
        Column name to write cell-type labels into.
    celltype_broad_key
        Column name to write broad cell-type labels into.
    unknown_label
        Label applied to clusters absent from the mapping.
    unknown_color
        Hex color for the unknown label.
    copy
        Return a copy rather than modifying in-place.

    Returns
    -------
    AnnData
        Modified adata (same object unless copy=True).
    """
    if copy:
        adata = adata.copy()

    # Load mapping if given as a path
    if not isinstance(mapping, dict):
        with open(mapping) as fh:
            mapping = json.load(fh)

    # Normalise keys to strings
    mapping = {str(k): v for k, v in mapping.items()}

    clusters = adata.obs[cluster_key].astype(str)

    # Build label and broad arrays
    celltype_arr = np.array(
        [mapping[c]["celltype"] if c in mapping else unknown_label for c in clusters],
        dtype=object,
    )
    celltype_broad_arr = np.array(
        [
            mapping[c].get("celltype_broad", mapping[c]["celltype"])
            if c in mapping
            else unknown_label
            for c in clusters
        ],
        dtype=object,
    )

    adata.obs[celltype_key] = pd.Categorical(celltype_arr)
    adata.obs[celltype_broad_key] = pd.Categorical(celltype_broad_arr)

    # Build color mapping for celltype_key
    categories = list(adata.obs[celltype_key].cat.categories)
    color_map: dict[str, str] = {}
    # Collect colors from the mapping dict first
    for _cluster_val, info in mapping.items():
        ct = info.get("celltype", unknown_label)
        if "color" in info:
            color_map[ct] = info["color"]

    # Assign palette colors to any category without an explicit color
    palette_iter = iter(_DEFAULT_PALETTE)
    for cat in categories:
        if cat == unknown_label:
            color_map.setdefault(cat, unknown_color)
        else:
            if cat not in color_map:
                try:
                    color_map[cat] = next(palette_iter)
                except StopIteration:
                    color_map[cat] = unknown_color

    colors_list = [color_map.get(cat, unknown_color) for cat in categories]
    colors_key = f"{celltype_key}_colors"
    # Only write if length differs (standard sc_tools convention)
    if colors_key not in adata.uns or len(adata.uns[colors_key]) != len(categories):
        adata.uns[colors_key] = colors_list

    return adata
