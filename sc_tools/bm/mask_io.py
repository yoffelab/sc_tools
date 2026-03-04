"""
Mask format adapters for segmentation benchmarking.

Loads segmentation masks from various tools into standardized labeled
numpy arrays (0=background, >0=cell_id).
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

__all__ = [
    "load_mask",
    "load_cellpose_mask",
    "load_stardist_mask",
    "load_cellprofiler_mask",
    "load_deepcell_mask",
    "load_tiff_mask",
]

logger = logging.getLogger(__name__)


def load_mask(path: str | Path, format: str = "auto") -> np.ndarray:
    """Load a segmentation mask, auto-detecting format from extension.

    Parameters
    ----------
    path
        Path to the mask file.
    format
        One of ``"auto"``, ``"cellpose"``, ``"stardist"``, ``"cellprofiler"``,
        ``"deepcell"``, ``"tiff"``. ``"auto"`` detects from file extension.

    Returns
    -------
    Labeled integer array (0=background, >0=cell_id).
    """
    path = Path(path)
    if format == "auto":
        format = _detect_format(path)

    loaders = {
        "cellpose": load_cellpose_mask,
        "stardist": load_stardist_mask,
        "cellprofiler": load_cellprofiler_mask,
        "deepcell": load_deepcell_mask,
        "tiff": load_tiff_mask,
    }
    if format not in loaders:
        raise ValueError(f"Unknown mask format: {format!r}. Choose from {list(loaders)}")
    return loaders[format](path)


def _detect_format(path: Path) -> str:
    """Detect mask format from file extension and naming conventions."""
    name = path.name.lower()
    suffix = path.suffix.lower()

    if name.endswith("_seg.npy") or suffix == ".npy":
        return "cellpose"
    if suffix in (".tif", ".tiff"):
        return "tiff"
    if suffix == ".npz":
        return "stardist"

    raise ValueError(
        f"Cannot auto-detect mask format for {path.name!r}. Specify format explicitly."
    )


def load_cellpose_mask(path: str | Path) -> np.ndarray:
    """Load a Cellpose segmentation mask from ``*_seg.npy``.

    Cellpose saves a dict with key ``"masks"`` containing the labeled array.
    Also handles plain labeled arrays saved directly.
    """
    path = Path(path)
    data = np.load(path, allow_pickle=True)

    if isinstance(data, np.ndarray) and data.ndim == 0:
        # dict saved via np.save
        item = data.item()
        if isinstance(item, dict) and "masks" in item:
            mask = np.asarray(item["masks"])
        else:
            raise ValueError(f"Cellpose file does not contain 'masks' key: {path}")
    elif isinstance(data, np.ndarray) and data.ndim >= 2:
        mask = data
    else:
        raise ValueError(f"Unexpected Cellpose file structure: {path}")

    return _validate_mask(mask, path)


def load_stardist_mask(path: str | Path) -> np.ndarray:
    """Load a StarDist label image from ``.tif``/``.tiff`` or ``.npz``.

    For ``.npz`` files, expects a ``"labels"`` key.
    """
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix == ".npz":
        data = np.load(path)
        if "labels" in data:
            mask = data["labels"]
        else:
            # Try first array
            keys = list(data.keys())
            if keys:
                mask = data[keys[0]]
            else:
                raise ValueError(f"Empty npz file: {path}")
    elif suffix in (".tif", ".tiff"):
        mask = _load_tiff(path)
    else:
        raise ValueError(f"Unsupported StarDist file extension: {suffix}")

    return _validate_mask(mask, path)


def load_cellprofiler_mask(path: str | Path) -> np.ndarray:
    """Load a CellProfiler label image from TIFF."""
    return load_tiff_mask(path)


def load_deepcell_mask(path: str | Path) -> np.ndarray:
    """Load a DeepCell/Mesmer output mask from TIFF or NPZ.

    DeepCell output is typically ``(1, H, W, 1)``; squeeze extra dims.
    """
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix == ".npz":
        data = np.load(path)
        keys = list(data.keys())
        mask = data[keys[0]] if keys else np.array([])
    elif suffix in (".tif", ".tiff"):
        mask = _load_tiff(path)
    else:
        mask = np.load(path)

    mask = np.squeeze(mask)
    return _validate_mask(mask, path)


def load_tiff_mask(path: str | Path) -> np.ndarray:
    """Load a labeled TIFF mask (generic loader)."""
    path = Path(path)
    mask = _load_tiff(path)
    return _validate_mask(mask, path)


def _load_tiff(path: Path) -> np.ndarray:
    """Load a TIFF file using tifffile."""
    try:
        import tifffile
    except ImportError as e:
        raise ImportError(
            "tifffile is required for TIFF mask loading. "
            "Install with: pip install sc-tools[benchmark]"
        ) from e
    return tifffile.imread(str(path))


def _validate_mask(mask: np.ndarray, path: Path) -> np.ndarray:
    """Validate and standardize a mask array."""
    if mask.ndim != 2:
        if mask.ndim > 2:
            mask = np.squeeze(mask)
        if mask.ndim != 2:
            raise ValueError(f"Mask must be 2D, got {mask.ndim}D array from {path}")

    # Ensure integer labels
    if not np.issubdtype(mask.dtype, np.integer):
        mask = mask.astype(np.int32)

    n_cells = len(np.unique(mask)) - (1 if 0 in mask else 0)
    logger.debug("Loaded mask from %s: shape=%s, n_cells=%d", path, mask.shape, n_cells)
    return mask
