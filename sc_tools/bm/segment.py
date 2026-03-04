"""
Thin wrappers for running Cellpose and StarDist on IMC multi-channel TIFFs.

Both functions accept a multi-channel intensity TIFF (C, H, W) and return
a labeled segmentation mask (H, W).
"""

from __future__ import annotations

import logging

import numpy as np

__all__ = ["run_cellpose", "run_stardist"]

logger = logging.getLogger(__name__)


def _extract_and_normalize(
    intensity: np.ndarray,
    channel_indices: list[int],
) -> np.ndarray:
    """Extract channels from a (C, H, W) array, average, and normalize to [0, 1]."""
    selected = intensity[channel_indices].astype(np.float32)
    combined = np.mean(selected, axis=0)
    vmin, vmax = combined.min(), combined.max()
    if vmax > vmin:
        combined = (combined - vmin) / (vmax - vmin)
    return combined


def run_cellpose(
    intensity_tiff: np.ndarray,
    nuclear_channels: list[int],
    membrane_channels: list[int] | None = None,
    model_type: str = "cyto2",
    diameter: float | None = None,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.0,
    gpu: bool = False,
) -> np.ndarray:
    """Run Cellpose segmentation on an IMC multi-channel TIFF.

    Parameters
    ----------
    intensity_tiff
        Multi-channel image, shape ``(C, H, W)``.
    nuclear_channels
        Indices of nuclear channels (e.g. ``[39, 40]`` for DNA1, DNA2).
    membrane_channels
        Indices of membrane channels (e.g. ``[14]`` for CD45). If None,
        runs nucleus-only (cyto2 still uses the nuclear channel as ch2).
    model_type
        Cellpose model type (default ``"cyto2"``).
    diameter
        Expected cell diameter in pixels. None = auto-estimate.
    flow_threshold
        Flow error threshold for Cellpose.
    cellprob_threshold
        Cell probability threshold for Cellpose.
    gpu
        Whether to use GPU.

    Returns
    -------
    Labeled segmentation mask, shape ``(H, W)``, dtype uint32.
    """
    try:
        from cellpose import models
    except ImportError as e:
        raise ImportError(
            "cellpose is required. Install with: pip install 'sc-tools[benchmark]'"
        ) from e

    nuclear = _extract_and_normalize(intensity_tiff, nuclear_channels)

    if membrane_channels is not None and len(membrane_channels) > 0:
        membrane = _extract_and_normalize(intensity_tiff, membrane_channels)
        # Cellpose expects shape (2, H, W): channel 0 = cytoplasm/membrane, channel 1 = nucleus
        img = np.stack([membrane, nuclear], axis=0)
        channels = [1, 2]  # [cytoplasm, nucleus]
    else:
        img = nuclear
        channels = [0, 0]  # grayscale

    # Support cellpose v3 (models.Cellpose) and v4+ (models.CellposeModel)
    try:
        model_cls = models.Cellpose
    except AttributeError:
        model_cls = models.CellposeModel
    model = model_cls(model_type=model_type, gpu=gpu)
    mask, _flows, _styles, _diams = model.eval(
        img,
        diameter=diameter,
        channels=channels,
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
    )

    logger.info("Cellpose: detected %d cells", len(np.unique(mask)) - 1)
    return mask.astype(np.uint32)


def run_stardist(
    intensity_tiff: np.ndarray,
    nuclear_channels: list[int],
    model_name: str = "2D_versatile_fluo",
    prob_thresh: float | None = None,
    nms_thresh: float | None = None,
    scale: float | None = None,
) -> np.ndarray:
    """Run StarDist segmentation on an IMC multi-channel TIFF.

    Parameters
    ----------
    intensity_tiff
        Multi-channel image, shape ``(C, H, W)``.
    nuclear_channels
        Indices of nuclear channels (e.g. ``[39, 40]`` for DNA1, DNA2).
    model_name
        StarDist pretrained model name (default ``"2D_versatile_fluo"``).
    prob_thresh
        Probability threshold. None = model default.
    nms_thresh
        Non-maximum suppression threshold. None = model default.
    scale
        Scale factor for the image (useful if pixel size differs from
        training data). None = no rescaling.

    Returns
    -------
    Labeled segmentation mask, shape ``(H, W)``, dtype uint32.
    """
    try:
        from stardist.models import StarDist2D
    except ImportError as e:
        raise ImportError(
            "stardist is required. Install with: pip install 'sc-tools[benchmark]'"
        ) from e

    nuclear = _extract_and_normalize(intensity_tiff, nuclear_channels)

    model = StarDist2D.from_pretrained(model_name)

    predict_kwargs = {}
    if prob_thresh is not None:
        predict_kwargs["prob_thresh"] = prob_thresh
    if nms_thresh is not None:
        predict_kwargs["nms_thresh"] = nms_thresh
    if scale is not None:
        predict_kwargs["scale"] = scale

    mask, _details = model.predict_instances(nuclear, **predict_kwargs)

    logger.info("StarDist: detected %d cells", len(np.unique(mask)) - 1)
    return mask.astype(np.uint32)
