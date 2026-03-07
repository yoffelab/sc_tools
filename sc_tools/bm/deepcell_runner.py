"""DeepCell/Mesmer wrapper for IMC segmentation benchmark.

Mesmer is a pre-trained deep learning model for nuclear and whole-cell
segmentation, part of the DeepCell ecosystem. It expects a 4D input
``(batch, H, W, 2)`` where channel 0 = nuclear, channel 1 = membrane.
"""

from __future__ import annotations

import logging

import numpy as np

__all__ = ["run_deepcell"]

logger = logging.getLogger(__name__)


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

    Parameters
    ----------
    image
        Either a probability map ``(H, W, C)`` from Ilastik (3 channels:
        background, nucleus, cytoplasm), or a multi-channel intensity TIFF
        ``(C, H, W)``. Detected automatically from shape.
    nuclear_channels
        For ``(C, H, W)`` input: indices of nuclear channels.
    membrane_channels
        For ``(C, H, W)`` input: indices of membrane channels.
    nuclear_idx
        For ``(H, W, C)`` probability maps: index of the nuclear channel.
    cytoplasm_idx
        For ``(H, W, C)`` probability maps: index of the cytoplasm channel.
    compartment
        ``"whole-cell"`` or ``"nuclear"``.
    image_mpp
        Microns per pixel (IMC default: 1.0).
    postprocess_kwargs
        Extra kwargs for ``app.predict()``.

    Returns
    -------
    Labeled segmentation mask, shape ``(H, W)``, dtype uint32.
    """
    try:
        from deepcell.applications import Mesmer
    except ImportError as e:
        raise ImportError(
            "deepcell is required. Install with: pip install 'sc-tools[benchmark-extended]'"
        ) from e

    nuclear, membrane = _prepare_inputs(
        image, nuclear_channels, membrane_channels, nuclear_idx, cytoplasm_idx
    )

    # Mesmer expects (batch, H, W, 2)
    combined = np.stack([nuclear, membrane], axis=-1)
    combined = combined[np.newaxis, ...]  # add batch dim

    app = Mesmer()

    predict_kwargs = {
        "image_mpp": image_mpp,
        "compartment": compartment,
    }
    if postprocess_kwargs:
        predict_kwargs["postprocess_kwargs"] = postprocess_kwargs

    mask = app.predict(combined, **predict_kwargs)

    # Output is (batch, H, W, 1) — squeeze to (H, W)
    mask = np.squeeze(mask)

    logger.info("DeepCell Mesmer: detected %d cells", len(np.unique(mask)) - 1)
    return mask.astype(np.uint32)


def _prepare_inputs(
    image: np.ndarray,
    nuclear_channels: list[int] | None,
    membrane_channels: list[int] | None,
    nuclear_idx: int,
    cytoplasm_idx: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Prepare nuclear and membrane channels for Mesmer.

    Returns (nuclear, membrane) both as (H, W) float32 arrays.
    """
    if image.ndim == 3 and image.shape[2] <= 4:
        # (H, W, C) probability map
        nuclear = _normalize(image[:, :, nuclear_idx])
        membrane = _normalize(image[:, :, cytoplasm_idx])
    elif image.ndim == 3 and nuclear_channels is not None:
        # (C, H, W) multi-channel intensity
        nuclear = _extract_and_normalize(image, nuclear_channels)
        if membrane_channels is not None and len(membrane_channels) > 0:
            membrane = _extract_and_normalize(image, membrane_channels)
        else:
            membrane = np.zeros_like(nuclear)
    elif image.ndim == 2:
        nuclear = _normalize(image)
        membrane = np.zeros_like(nuclear)
    else:
        raise ValueError(
            f"Unexpected image shape {image.shape}. Expected (H, W, C) probability map "
            f"or (C, H, W) intensity TIFF."
        )

    return nuclear, membrane


def _normalize(img: np.ndarray) -> np.ndarray:
    """Normalize to [0, 1] float32."""
    img = img.astype(np.float32)
    vmin, vmax = img.min(), img.max()
    if vmax > vmin:
        img = (img - vmin) / (vmax - vmin)
    return img


def _extract_and_normalize(
    intensity: np.ndarray,
    channel_indices: list[int],
) -> np.ndarray:
    """Extract channels from (C, H, W), average, normalize."""
    selected = intensity[channel_indices].astype(np.float32)
    combined = np.mean(selected, axis=0)
    return _normalize(combined)
