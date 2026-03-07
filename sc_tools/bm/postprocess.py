"""Shared post-processing utilities for segmentation masks.

Converts semantic segmentation maps to instance masks and applies
filtering/cleanup operations.
"""

from __future__ import annotations

import logging

import numpy as np

__all__ = [
    "semantic_to_instance",
    "filter_masks_by_area",
    "fill_holes",
]

logger = logging.getLogger(__name__)


def semantic_to_instance(
    semantic_mask: np.ndarray,
    foreground_class: int = 1,
    min_distance: int = 5,
    compactness: float = 0.0,
) -> np.ndarray:
    """Convert a semantic segmentation map to an instance segmentation mask.

    Uses distance transform + watershed to split connected foreground
    regions into individual cells.

    Parameters
    ----------
    semantic_mask
        2D array with class labels (0=background, foreground_class=cells).
    foreground_class
        Which class label represents foreground/cells.
    min_distance
        Minimum distance between cell centers for peak detection.
    compactness
        Compactness parameter for watershed (0 = standard).

    Returns
    -------
    Labeled instance mask (H, W), dtype int32.
    """
    from scipy import ndimage
    from skimage.feature import peak_local_max
    from skimage.segmentation import watershed

    binary = (semantic_mask == foreground_class).astype(np.uint8)

    if np.sum(binary) == 0:
        return np.zeros_like(semantic_mask, dtype=np.int32)

    # Distance transform
    distance = ndimage.distance_transform_edt(binary)

    # Find peaks (cell centers)
    coords = peak_local_max(
        distance,
        min_distance=min_distance,
        labels=binary,
    )

    if len(coords) == 0:
        # Single connected component
        labels, _ = ndimage.label(binary)
        return labels.astype(np.int32)

    # Create markers from peaks
    markers = np.zeros_like(semantic_mask, dtype=np.int32)
    for i, (r, c) in enumerate(coords, 1):
        markers[r, c] = i

    # Watershed
    labels = watershed(
        -distance,
        markers=markers,
        mask=binary,
        compactness=compactness,
    )

    logger.debug(
        "semantic_to_instance: %d instances from semantic mask", len(np.unique(labels)) - 1
    )
    return labels.astype(np.int32)


def filter_masks_by_area(
    mask: np.ndarray,
    min_area: int = 10,
    max_area: int | None = None,
) -> np.ndarray:
    """Remove cells outside area bounds.

    Parameters
    ----------
    mask
        Labeled instance mask.
    min_area
        Minimum cell area in pixels.
    max_area
        Maximum cell area in pixels. None = no upper limit.

    Returns
    -------
    Filtered mask with cells relabeled contiguously.
    """
    from skimage.measure import regionprops

    props = regionprops(mask)
    keep_labels = set()

    for prop in props:
        if prop.area < min_area:
            continue
        if max_area is not None and prop.area > max_area:
            continue
        keep_labels.add(prop.label)

    if len(keep_labels) == len(props):
        return mask

    # Build filtered mask with contiguous labels
    filtered = np.zeros_like(mask)
    new_label = 1
    for old_label in sorted(keep_labels):
        filtered[mask == old_label] = new_label
        new_label += 1

    n_removed = len(props) - len(keep_labels)
    logger.debug("filter_masks_by_area: removed %d cells, kept %d", n_removed, len(keep_labels))
    return filtered


def fill_holes(mask: np.ndarray, max_hole_area: int = 50) -> np.ndarray:
    """Fill small holes within cell regions.

    Parameters
    ----------
    mask
        Labeled instance mask.
    max_hole_area
        Maximum hole area to fill.

    Returns
    -------
    Mask with small holes filled.
    """
    from scipy import ndimage

    result = mask.copy()
    labels = np.unique(mask)
    labels = labels[labels > 0]

    for lbl in labels:
        cell_mask = mask == lbl
        filled = ndimage.binary_fill_holes(cell_mask)
        # Only fill if the new area is not too much larger
        hole_area = np.sum(filled & ~cell_mask)
        if hole_area <= max_hole_area:
            result[filled & ~cell_mask] = lbl

    return result
