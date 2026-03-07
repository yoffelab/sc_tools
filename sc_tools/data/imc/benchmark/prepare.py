"""ROI validation, probability map generation, channel extraction, and normalization.

Handles the critical IMC intensity normalization: raw pixel values are ion counts
(typically 0-100+ range, float32), NOT standard 0-255 uint8 images.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

__all__ = [
    "validate_roi_files",
    "generate_probability_map",
    "extract_dna_channels",
    "extract_channel_by_name",
    "normalize_imc_intensity",
]

logger = logging.getLogger(__name__)

# Common DNA channel names in IMC panels
DNA_CHANNEL_NAMES = [
    "DNA1",
    "DNA2",
    "DNA3",
    "Ir191",
    "Ir193",  # Iridium intercalators
    "191Ir",
    "193Ir",
    "Histone_H3",
    "HistoneH3",
]


def validate_roi_files(
    tiff_path: str | Path,
    channel_csv_path: str | Path | None = None,
) -> dict[str, bool]:
    """Validate that an ROI has required files and correct format.

    Parameters
    ----------
    tiff_path
        Path to the ``*_full.tiff`` file.
    channel_csv_path
        Optional path to the ``*_full.csv`` channel index file.

    Returns
    -------
    Dict with keys: ``tiff_exists``, ``tiff_readable``, ``has_channels``,
    ``csv_exists``, ``channels_match``, ``valid`` (all checks pass).
    """
    tiff_path = Path(tiff_path)
    result = {
        "tiff_exists": tiff_path.is_file(),
        "tiff_readable": False,
        "has_channels": False,
        "csv_exists": False,
        "channels_match": False,
        "valid": False,
    }

    if not result["tiff_exists"]:
        return result

    try:
        import tifffile

        img = tifffile.imread(str(tiff_path))
        result["tiff_readable"] = True
        result["has_channels"] = img.ndim == 3 and img.shape[0] > 1
        n_channels_tiff = img.shape[0] if img.ndim == 3 else 1
    except Exception as e:
        logger.warning("Cannot read TIFF %s: %s", tiff_path, e)
        return result

    if channel_csv_path is not None:
        channel_csv_path = Path(channel_csv_path)
        result["csv_exists"] = channel_csv_path.is_file()
        if result["csv_exists"]:
            try:
                import pandas as pd

                csv_df = pd.read_csv(channel_csv_path, header=None)
                n_channels_csv = len(csv_df)
                result["channels_match"] = n_channels_csv == n_channels_tiff
            except Exception as e:
                logger.warning("Cannot read channel CSV %s: %s", channel_csv_path, e)
    else:
        # No CSV to validate against
        result["csv_exists"] = True
        result["channels_match"] = True

    result["valid"] = all(
        [
            result["tiff_exists"],
            result["tiff_readable"],
            result["has_channels"],
            result["channels_match"],
        ]
    )
    return result


def normalize_imc_intensity(
    image: np.ndarray,
    method: str = "percentile",
    low_pct: float = 1.0,
    high_pct: float = 99.5,
    cofactor: float = 5.0,
) -> np.ndarray:
    """Normalize raw IMC ion count intensities.

    Parameters
    ----------
    image
        Raw IMC image, shape ``(H, W)`` or ``(C, H, W)``.
    method
        Normalization method:
        - ``"percentile"``: clip to [low_pct, high_pct] percentile, scale to [0, 1]
        - ``"arcsinh"``: arcsinh(x / cofactor), standard for IMC
        - ``"zscore"``: per-channel z-score (mean=0, std=1)
        - ``"uint8"``: percentile clip + scale to [0, 255] uint8
    low_pct
        Lower percentile for clipping (percentile/uint8 methods).
    high_pct
        Upper percentile for clipping (percentile/uint8 methods).
    cofactor
        Cofactor for arcsinh transform.

    Returns
    -------
    Normalized image (same shape as input).
    """
    image = np.asarray(image, dtype=np.float32)

    if method == "percentile":
        return _normalize_percentile(image, low_pct, high_pct)
    elif method == "arcsinh":
        return np.arcsinh(image / cofactor)
    elif method == "zscore":
        return _normalize_zscore(image)
    elif method == "uint8":
        return _normalize_uint8(image, low_pct, high_pct)
    else:
        raise ValueError(f"Unknown normalization method: {method!r}")


def _normalize_percentile(image: np.ndarray, low_pct: float, high_pct: float) -> np.ndarray:
    """Percentile-based normalization to [0, 1]."""
    if image.ndim == 2:
        vmin = np.percentile(image, low_pct)
        vmax = np.percentile(image, high_pct)
        if vmax > vmin:
            return np.clip((image - vmin) / (vmax - vmin), 0, 1)
        return np.zeros_like(image)

    # Per-channel normalization for (C, H, W)
    result = np.zeros_like(image)
    for c in range(image.shape[0]):
        result[c] = _normalize_percentile(image[c], low_pct, high_pct)
    return result


def _normalize_zscore(image: np.ndarray) -> np.ndarray:
    """Per-channel z-score normalization."""
    if image.ndim == 2:
        mu = image.mean()
        std = image.std()
        return (image - mu) / max(std, 1e-8)

    result = np.zeros_like(image)
    for c in range(image.shape[0]):
        result[c] = _normalize_zscore(image[c])
    return result


def _normalize_uint8(image: np.ndarray, low_pct: float, high_pct: float) -> np.ndarray:
    """Percentile clip + scale to [0, 255] uint8."""
    normed = _normalize_percentile(image, low_pct, high_pct)
    return (normed * 255).astype(np.uint8)


def extract_dna_channels(
    tiff_path: str | Path,
    channel_csv_path: str | Path | None = None,
    panel_mapper=None,
) -> np.ndarray:
    """Extract DNA channels from a multi-channel IMC TIFF.

    Parameters
    ----------
    tiff_path
        Path to ``*_full.tiff``.
    channel_csv_path
        Path to ``*_full.csv`` channel index.
    panel_mapper
        Optional ``IMCPanelMapper`` instance. If None, one is created from
        the channel CSV.

    Returns
    -------
    2D array (H, W) — averaged DNA signal from all found DNA channels.
    """
    import tifffile

    image = tifffile.imread(str(tiff_path))
    if image.ndim != 3:
        raise ValueError(f"Expected 3D TIFF (C, H, W), got shape {image.shape}")

    if panel_mapper is None:
        panel_mapper = _make_panel_mapper(channel_csv_path)

    # Find DNA channel indices
    dna_indices = []
    for name in DNA_CHANNEL_NAMES:
        idx = panel_mapper.resolve(name)
        if idx is not None and idx < image.shape[0]:
            dna_indices.append(idx)

    if not dna_indices:
        logger.warning("No DNA channels found. Using first channel as fallback.")
        dna_indices = [0]

    # Average DNA channels
    dna_stack = image[dna_indices].astype(np.float32)
    return np.mean(dna_stack, axis=0)


def extract_channel_by_name(
    tiff_path: str | Path,
    channel_name: str,
    channel_csv_path: str | Path | None = None,
    panel_mapper=None,
) -> np.ndarray | None:
    """Extract a single channel from a multi-channel IMC TIFF by name.

    Parameters
    ----------
    tiff_path
        Path to ``*_full.tiff``.
    channel_name
        Protein name, isotope tag, or full ``MarkerName(IsotopeTag)`` string.
    channel_csv_path
        Path to ``*_full.csv`` channel index.
    panel_mapper
        Optional ``IMCPanelMapper`` instance.

    Returns
    -------
    2D array (H, W) or None if channel not found.
    """
    import tifffile

    image = tifffile.imread(str(tiff_path))
    if image.ndim != 3:
        raise ValueError(f"Expected 3D TIFF (C, H, W), got shape {image.shape}")

    if panel_mapper is None:
        panel_mapper = _make_panel_mapper(channel_csv_path)

    idx = panel_mapper.resolve(channel_name)
    if idx is None or idx >= image.shape[0]:
        logger.warning("Channel %r not found or out of range", channel_name)
        return None

    return image[idx].astype(np.float32)


def generate_probability_map(
    dna_image: np.ndarray,
    method: str = "gaussian",
    sigma: float = 2.0,
) -> np.ndarray:
    """Generate a probability map from a DNA channel image (no Ilastik needed).

    Parameters
    ----------
    dna_image
        2D DNA channel image (H, W).
    method
        Generation method:
        - ``"gaussian"``: Gaussian smoothing + Otsu threshold → 3-channel prob map
        - ``"otsu"``: Simple Otsu threshold → binary prob map
        - ``"multiscale"``: Multi-scale Gaussian → probability estimate

    Returns
    -------
    Probability map (H, W, 3): [background, nucleus, cytoplasm].
    """
    from scipy.ndimage import gaussian_filter

    dna = dna_image.astype(np.float32)

    # Normalize to [0, 1]
    vmin, vmax = np.percentile(dna, [1, 99.5])
    if vmax > vmin:
        dna = np.clip((dna - vmin) / (vmax - vmin), 0, 1)

    if method == "gaussian":
        smoothed = gaussian_filter(dna, sigma=sigma)

        # Otsu-like threshold
        try:
            from skimage.filters import threshold_otsu

            thresh = threshold_otsu(smoothed[smoothed > 0]) if np.any(smoothed > 0) else 0.5
        except (ImportError, ValueError):
            thresh = 0.5

        nuclear_prob = np.clip(smoothed / max(thresh * 2, 1e-8), 0, 1)
        bg_prob = 1.0 - nuclear_prob

        # Estimate cytoplasm as ring around nuclei
        from scipy.ndimage import binary_dilation

        nuclear_mask = smoothed > thresh
        dilated = binary_dilation(nuclear_mask, iterations=3)
        cyto_mask = dilated & ~nuclear_mask
        cyto_prob = np.where(cyto_mask, 0.7, 0.1)

        # Normalize probabilities
        total = bg_prob + nuclear_prob + cyto_prob
        prob_map = np.stack([bg_prob / total, nuclear_prob / total, cyto_prob / total], axis=-1)

    elif method == "otsu":
        try:
            from skimage.filters import threshold_otsu

            thresh = threshold_otsu(dna[dna > 0]) if np.any(dna > 0) else 0.5
        except (ImportError, ValueError):
            thresh = 0.5

        nuclear = (dna > thresh).astype(np.float32)
        bg = 1.0 - nuclear
        cyto = np.zeros_like(dna)
        prob_map = np.stack([bg, nuclear, cyto], axis=-1)

    elif method == "multiscale":
        # Multi-scale Gaussian smoothing for better probability estimation
        sigmas = [1.0, 2.0, 4.0]
        smoothed_stack = np.stack([gaussian_filter(dna, s) for s in sigmas])
        nuclear_prob = np.clip(np.mean(smoothed_stack, axis=0), 0, 1)
        bg_prob = 1.0 - nuclear_prob
        cyto_prob = np.clip(nuclear_prob * 0.3, 0, 0.5)

        total = bg_prob + nuclear_prob + cyto_prob
        prob_map = np.stack([bg_prob / total, nuclear_prob / total, cyto_prob / total], axis=-1)
    else:
        raise ValueError(f"Unknown probability map method: {method!r}")

    return prob_map.astype(np.float32)


def _make_panel_mapper(channel_csv_path: str | Path | None = None):
    """Create an IMCPanelMapper from a channel CSV."""
    from sc_tools.ingest.imc import IMCPanelMapper

    mapper = IMCPanelMapper()
    if channel_csv_path is not None:
        mapper.from_full_csv(str(channel_csv_path))
    return mapper
