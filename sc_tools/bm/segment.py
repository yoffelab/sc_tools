"""
Thin wrappers for running Cellpose and StarDist on IMC data.

Both functions accept either:
- A probability map ``(H, W, C)`` from Ilastik pixel classification (typical
  IMC pipeline: ``*_Probabilities.tiff`` with 3 channels: background, nucleus,
  cytoplasm).
- A multi-channel intensity TIFF ``(C, H, W)`` with explicit channel indices.

Both return a labeled segmentation mask ``(H, W)``.
"""

from __future__ import annotations

import logging

import numpy as np

__all__ = ["run_cellpose", "run_stardist", "run_deepcell", "run_all_strategy1"]

logger = logging.getLogger(__name__)


def _normalize(img: np.ndarray) -> np.ndarray:
    """Normalize a 2D image to [0, 1] float32."""
    img = img.astype(np.float32)
    vmin, vmax = img.min(), img.max()
    if vmax > vmin:
        img = (img - vmin) / (vmax - vmin)
    return img


def _extract_and_normalize(
    intensity: np.ndarray,
    channel_indices: list[int],
) -> np.ndarray:
    """Extract channels from a (C, H, W) array, average, and normalize to [0, 1]."""
    selected = intensity[channel_indices].astype(np.float32)
    combined = np.mean(selected, axis=0)
    return _normalize(combined)


def run_cellpose(
    image: np.ndarray,
    nuclear_channels: list[int] | None = None,
    membrane_channels: list[int] | None = None,
    nuclear_idx: int = 1,
    cytoplasm_idx: int = 2,
    model_type: str = "cyto2",
    diameter: float | None = None,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.0,
    gpu: bool = False,
) -> np.ndarray:
    """Run Cellpose segmentation.

    Parameters
    ----------
    image
        Either a probability map ``(H, W, C)`` from Ilastik (e.g. 3 channels:
        background, nucleus, cytoplasm), or a multi-channel intensity TIFF
        ``(C, H, W)``. Detected automatically from shape.
    nuclear_channels
        For ``(C, H, W)`` input: indices of nuclear channels. Ignored for
        ``(H, W, C)`` probability maps.
    membrane_channels
        For ``(C, H, W)`` input: indices of membrane channels. Ignored for
        ``(H, W, C)`` probability maps.
    nuclear_idx
        For ``(H, W, C)`` probability maps: index of the nuclear channel
        (default 1).
    cytoplasm_idx
        For ``(H, W, C)`` probability maps: index of the cytoplasm channel
        (default 2).
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

    # Detect input format
    if image.ndim == 3 and image.shape[2] <= 4:
        # (H, W, C) probability map
        nuclear = _normalize(image[:, :, nuclear_idx])
        cytoplasm = _normalize(image[:, :, cytoplasm_idx])
        img = np.stack([cytoplasm, nuclear], axis=0)
        channels = [1, 2]
    elif image.ndim == 3 and nuclear_channels is not None:
        # (C, H, W) multi-channel intensity
        nuclear = _extract_and_normalize(image, nuclear_channels)
        if membrane_channels is not None and len(membrane_channels) > 0:
            membrane = _extract_and_normalize(image, membrane_channels)
            img = np.stack([membrane, nuclear], axis=0)
            channels = [1, 2]
        else:
            img = nuclear
            channels = [0, 0]
    elif image.ndim == 2:
        img = _normalize(image)
        channels = [0, 0]
    else:
        raise ValueError(
            f"Unexpected image shape {image.shape}. Expected (H, W, C) probability map "
            f"or (C, H, W) intensity TIFF."
        )

    # Support cellpose v3 (models.Cellpose) and v4+ (models.CellposeModel)
    try:
        model_cls = models.Cellpose
    except AttributeError:
        model_cls = models.CellposeModel
    model = model_cls(model_type=model_type, gpu=gpu)
    result = model.eval(
        img,
        diameter=diameter,
        channels=channels,
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
    )
    # Cellpose v3 returns 4 values, v4+ returns 3 (no diams)
    mask = result[0]

    logger.info("Cellpose: detected %d cells", len(np.unique(mask)) - 1)
    return mask.astype(np.uint32)


def run_stardist(
    image: np.ndarray,
    nuclear_channels: list[int] | None = None,
    nuclear_idx: int = 1,
    model_name: str = "2D_versatile_fluo",
    prob_thresh: float | None = None,
    nms_thresh: float | None = None,
    scale: float | None = None,
) -> np.ndarray:
    """Run StarDist segmentation.

    Parameters
    ----------
    image
        Either a probability map ``(H, W, C)`` from Ilastik (nuclear channel
        at ``nuclear_idx``), or a multi-channel intensity TIFF ``(C, H, W)``
        (nuclear channels at ``nuclear_channels``), or a 2D nuclear image
        ``(H, W)``.
    nuclear_channels
        For ``(C, H, W)`` input: indices of nuclear channels.
    nuclear_idx
        For ``(H, W, C)`` probability maps: index of the nuclear channel
        (default 1).
    model_name
        StarDist pretrained model name (default ``"2D_versatile_fluo"``).
    prob_thresh
        Probability threshold. None = model default.
    nms_thresh
        Non-maximum suppression threshold. None = model default.
    scale
        Scale factor for the image. None = no rescaling.

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

    # Detect input format
    if image.ndim == 3 and image.shape[2] <= 4:
        # (H, W, C) probability map
        nuclear = _normalize(image[:, :, nuclear_idx])
    elif image.ndim == 3 and nuclear_channels is not None:
        # (C, H, W) multi-channel intensity
        nuclear = _extract_and_normalize(image, nuclear_channels)
    elif image.ndim == 2:
        nuclear = _normalize(image)
    else:
        raise ValueError(
            f"Unexpected image shape {image.shape}. Expected (H, W, C) probability map, "
            f"(C, H, W) intensity TIFF, or (H, W) nuclear image."
        )

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


def run_deepcell(
    image: np.ndarray,
    nuclear_channels: list[int] | None = None,
    membrane_channels: list[int] | None = None,
    nuclear_idx: int = 1,
    cytoplasm_idx: int = 2,
    compartment: str = "whole-cell",
    image_mpp: float = 1.0,
    postprocess_kwargs: dict | None = None,
) -> np.ndarray:
    """Run DeepCell Mesmer segmentation.

    Thin re-export of ``sc_tools.bm.deepcell_runner.run_deepcell``.
    See that module for full documentation.

    Returns
    -------
    Labeled segmentation mask, shape ``(H, W)``, dtype uint32.
    """
    from sc_tools.bm.deepcell_runner import run_deepcell as _run_deepcell

    return _run_deepcell(
        image,
        nuclear_channels=nuclear_channels,
        membrane_channels=membrane_channels,
        nuclear_idx=nuclear_idx,
        cytoplasm_idx=cytoplasm_idx,
        compartment=compartment,
        image_mpp=image_mpp,
        postprocess_kwargs=postprocess_kwargs,
    )


def run_all_strategy1(
    image: np.ndarray,
    methods: list[str] | None = None,
    nuclear_channels: list[int] | None = None,
    membrane_channels: list[int] | None = None,
    nuclear_idx: int = 1,
    cytoplasm_idx: int = 2,
    gpu: bool = False,
) -> dict[str, np.ndarray]:
    """Run all Strategy 1 methods on a single image (Ilastik prob map or intensity TIFF).

    Parameters
    ----------
    image
        Probability map ``(H, W, C)`` or intensity TIFF ``(C, H, W)``.
    methods
        Method names to run. Default: all available.
    nuclear_channels
        Nuclear channel indices for ``(C, H, W)`` input.
    membrane_channels
        Membrane channel indices for ``(C, H, W)`` input.
    nuclear_idx
        Nuclear channel index for ``(H, W, C)`` input.
    cytoplasm_idx
        Cytoplasm channel index for ``(H, W, C)`` input.
    gpu
        Whether to use GPU.

    Returns
    -------
    Dict mapping method name to labeled mask.
    """
    if methods is None:
        methods = ["cellpose_cyto2", "cellpose_cyto3", "cellpose_nuclei", "stardist", "deepcell"]

    results = {}
    common = {
        "nuclear_channels": nuclear_channels,
        "membrane_channels": membrane_channels,
        "nuclear_idx": nuclear_idx,
        "cytoplasm_idx": cytoplasm_idx,
    }

    for method in methods:
        try:
            if method.startswith("cellpose"):
                model_type = method.replace("cellpose_", "") if "_" in method else "cyto2"
                mask = run_cellpose(image, model_type=model_type, gpu=gpu, **common)
            elif method == "stardist":
                mask = run_stardist(
                    image, **{k: v for k, v in common.items() if "membrane" not in k}
                )
            elif method == "deepcell":
                mask = run_deepcell(image, **common)
            else:
                logger.warning("Unknown Strategy 1 method: %s", method)
                continue
            results[f"s1_{method}"] = mask
            logger.info("Strategy 1 / %s: %d cells", method, len(np.unique(mask)) - 1)
        except Exception as e:
            logger.error("Strategy 1 / %s failed: %s", method, e)

    return results
