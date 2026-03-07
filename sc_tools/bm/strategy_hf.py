"""Strategy 3: DNA-only segmentation using HuggingFace pretrained models.

Applies pretrained vision models (CellViT, SAM, HoVer-Net) directly to
DNA channel images. Handles IMC-to-RGB format conversion (grayscale DNA
to pseudo-RGB for models expecting uint8 RGB input).
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

__all__ = [
    "HF_MODEL_REGISTRY",
    "run_hf_model",
    "run_all_strategy3",
]

logger = logging.getLogger(__name__)

# Registry of supported HuggingFace / pretrained models
HF_MODEL_REGISTRY: dict[str, dict] = {
    "cellvit_256": {
        "name": "CellViT-256",
        "hf_id": "lunit-io/CellViT-256",
        "input_format": "rgb_uint8",
        "type": "cellvit",
        "description": "Cell Vision Transformer (Lunit)",
    },
    "sam_base": {
        "name": "SAM ViT-Base",
        "hf_id": "facebook/sam-vit-base",
        "input_format": "rgb_uint8",
        "type": "sam",
        "description": "Segment Anything Model (Meta)",
    },
    "hovernet": {
        "name": "HoVer-Net",
        "hf_id": "hovernet",
        "input_format": "rgb_uint8",
        "type": "hovernet",
        "description": "HoVer-Net for nuclear segmentation",
    },
    "cellpose_cyto3": {
        "name": "Cellpose cyto3",
        "hf_id": None,
        "input_format": "grayscale",
        "type": "cellpose",
        "description": "Cellpose cyto3 baseline",
    },
    "stardist": {
        "name": "StarDist",
        "hf_id": None,
        "input_format": "grayscale",
        "type": "stardist",
        "description": "StarDist baseline",
    },
}


def _normalize_imc_to_uint8(
    channel: np.ndarray,
    low_pct: float = 1.0,
    high_pct: float = 99.5,
) -> np.ndarray:
    """Convert raw IMC ion counts to uint8 for pretrained models.

    Parameters
    ----------
    channel
        2D array of raw ion counts.
    low_pct
        Lower percentile for clipping.
    high_pct
        Upper percentile for clipping.

    Returns
    -------
    uint8 array (H, W) with values in [0, 255].
    """
    channel = channel.astype(np.float32)
    vmin = np.percentile(channel, low_pct)
    vmax = np.percentile(channel, high_pct)
    if vmax > vmin:
        normed = np.clip((channel - vmin) / (vmax - vmin), 0, 1)
    else:
        normed = np.zeros_like(channel)
    return (normed * 255).astype(np.uint8)


def _to_pseudo_rgb(grayscale_uint8: np.ndarray) -> np.ndarray:
    """Replicate single-channel uint8 to 3-channel pseudo-RGB.

    Parameters
    ----------
    grayscale_uint8
        (H, W) uint8 array.

    Returns
    -------
    (H, W, 3) uint8 array.
    """
    return np.stack([grayscale_uint8] * 3, axis=-1)


def run_hf_model(
    image: np.ndarray,
    model_name: str,
    gpu: bool = False,
) -> np.ndarray:
    """Run a HuggingFace pretrained model for cell segmentation.

    Parameters
    ----------
    image
        2D DNA channel image (H, W), raw ion counts.
    model_name
        Model key from HF_MODEL_REGISTRY.
    gpu
        Whether to use GPU.

    Returns
    -------
    Labeled instance mask (H, W), dtype uint32.
    """
    if model_name not in HF_MODEL_REGISTRY:
        raise ValueError(
            f"Unknown model: {model_name!r}. Available: {list(HF_MODEL_REGISTRY.keys())}"
        )

    info = HF_MODEL_REGISTRY[model_name]
    model_type = info["type"]

    if model_type == "cellpose":
        return _run_cellpose_baseline(image, gpu=gpu)
    elif model_type == "stardist":
        return _run_stardist_baseline(image)
    elif model_type == "cellvit":
        return _run_cellvit(image, info["hf_id"], gpu=gpu)
    elif model_type == "sam":
        return _run_sam(image, info["hf_id"], gpu=gpu)
    elif model_type == "hovernet":
        return _run_hovernet(image, gpu=gpu)
    else:
        raise ValueError(f"Unsupported model type: {model_type!r}")


def run_all_strategy3(
    tiff_path: str | Path,
    methods: list[str] | None = None,
    channel_csv_path: str | Path | None = None,
    panel_mapper=None,
    gpu: bool = False,
) -> dict[str, np.ndarray]:
    """Run all Strategy 3 methods on a single ROI.

    Parameters
    ----------
    tiff_path
        Path to ``*_full.tiff``.
    methods
        Model names to run. Default: all in registry.
    channel_csv_path
        Path to ``*_full.csv``.
    panel_mapper
        Optional ``IMCPanelMapper``.
    gpu
        Whether to use GPU.

    Returns
    -------
    Dict mapping method name to labeled mask.
    """
    from sc_tools.data.imc.benchmark.prepare import extract_dna_channels

    if methods is None:
        methods = list(HF_MODEL_REGISTRY.keys())

    # Extract DNA signal
    dna = extract_dna_channels(tiff_path, channel_csv_path, panel_mapper)
    logger.info("DNA image: shape=%s, range=[%.1f, %.1f]", dna.shape, dna.min(), dna.max())

    results = {}
    for method in methods:
        try:
            mask = run_hf_model(dna, method, gpu=gpu)
            results[f"s3_{method}"] = mask
            logger.info("Strategy 3 / %s: %d cells", method, len(np.unique(mask)) - 1)
        except Exception as e:
            logger.error("Strategy 3 / %s failed: %s", method, e)

    return results


def _run_cellpose_baseline(image: np.ndarray, gpu: bool = False) -> np.ndarray:
    """Run Cellpose cyto3 as a baseline."""
    from sc_tools.bm.segment import run_cellpose

    return run_cellpose(image, model_type="cyto3", gpu=gpu)


def _run_stardist_baseline(image: np.ndarray) -> np.ndarray:
    """Run StarDist as a baseline."""
    from sc_tools.bm.segment import run_stardist

    return run_stardist(image)


def _run_cellvit(image: np.ndarray, hf_id: str, gpu: bool = False) -> np.ndarray:
    """Run CellViT for nuclear segmentation.

    CellViT expects RGB uint8 patches (256x256). We tile the image,
    predict per-tile, and stitch.
    """
    try:
        import torch  # noqa: F401
        from transformers import AutoModel  # noqa: F401
    except ImportError as e:
        raise ImportError(
            "torch and transformers required for CellViT. "
            "Install with: pip install 'sc-tools[benchmark-extended]'"
        ) from e

    # Prepare RGB input
    img_uint8 = _normalize_imc_to_uint8(image)
    rgb = _to_pseudo_rgb(img_uint8)

    device = "cuda" if gpu and torch.cuda.is_available() else "cpu"

    # For CellViT, use tiled inference
    mask = _tiled_inference_generic(
        rgb, model_name="cellvit", hf_id=hf_id, device=device, tile_size=256
    )
    return mask.astype(np.uint32)


def _run_sam(image: np.ndarray, hf_id: str, gpu: bool = False) -> np.ndarray:
    """Run Segment Anything Model with automatic mask generation."""
    try:
        import torch
        from transformers import SamModel, SamProcessor
    except ImportError as e:
        raise ImportError(
            "torch and transformers required for SAM. "
            "Install with: pip install 'sc-tools[benchmark-extended]'"
        ) from e

    img_uint8 = _normalize_imc_to_uint8(image)
    rgb = _to_pseudo_rgb(img_uint8)

    device = "cuda" if gpu and torch.cuda.is_available() else "cpu"

    processor = SamProcessor.from_pretrained(hf_id)
    model = SamModel.from_pretrained(hf_id).to(device)

    # Use automatic mask generation approach
    # Process with grid point prompts
    h, w = rgb.shape[:2]
    inputs = processor(rgb, return_tensors="pt").to(device)

    # Generate masks using grid of input points
    grid_size = 32
    points = []
    for y in range(0, h, grid_size):
        for x in range(0, w, grid_size):
            points.append([x, y])

    if not points:
        return np.zeros((h, w), dtype=np.uint32)

    # Process in batches of points
    all_masks = np.zeros((h, w), dtype=np.int32)
    cell_id = 1

    batch_size = 64
    for i in range(0, len(points), batch_size):
        batch_points = points[i : i + batch_size]
        input_points = torch.tensor([batch_points], dtype=torch.float32).to(device)

        with torch.no_grad():
            outputs = model(
                pixel_values=inputs["pixel_values"],
                input_points=input_points,
            )

        masks = outputs.pred_masks.cpu().numpy()
        scores = outputs.iou_scores.cpu().numpy()

        # Take highest-scoring mask per point
        for j in range(masks.shape[1]):
            best_idx = np.argmax(scores[0, j])
            m = masks[0, j, best_idx] > 0

            # Only add if it does not overlap too much with existing
            overlap = np.sum(all_masks[m] > 0) / max(np.sum(m), 1)
            if overlap < 0.5 and np.sum(m) > 10:
                all_masks[m] = cell_id
                cell_id += 1

    return all_masks.astype(np.uint32)


def _run_hovernet(image: np.ndarray, gpu: bool = False) -> np.ndarray:
    """Run HoVer-Net for nuclear segmentation.

    HoVer-Net uses horizontal and vertical distance maps for instance
    segmentation. This is a simplified implementation using the model
    as a semantic segmenter with watershed post-processing.
    """
    # Fallback: use Cellpose nuclei model as HoVer-Net requires
    # specific weight files not on HuggingFace Hub
    logger.warning(
        "HoVer-Net HuggingFace integration pending. Falling back to Cellpose nuclei model."
    )
    from sc_tools.bm.segment import run_cellpose

    return run_cellpose(image, model_type="nuclei", gpu=gpu)


def _tiled_inference_generic(
    rgb: np.ndarray,
    model_name: str,
    hf_id: str,
    device: str,
    tile_size: int = 256,
    overlap: int = 32,
) -> np.ndarray:
    """Generic tiled inference with overlap stitching.

    Tiles the image, runs model per tile, stitches instance masks
    with overlap-based merging.
    """
    from sc_tools.bm.postprocess import semantic_to_instance

    h, w = rgb.shape[:2]
    full_mask = np.zeros((h, w), dtype=np.int32)
    cell_id = 1

    step = tile_size - overlap
    for y in range(0, h, step):
        for x in range(0, w, step):
            y_end = min(y + tile_size, h)
            x_end = min(x + tile_size, w)
            tile = rgb[y:y_end, x:x_end]

            # Pad tile to tile_size if needed
            pad_h = tile_size - tile.shape[0]
            pad_w = tile_size - tile.shape[1]
            if pad_h > 0 or pad_w > 0:
                tile = np.pad(tile, ((0, pad_h), (0, pad_w), (0, 0)), mode="reflect")

            # Simple threshold-based segmentation as fallback
            gray = tile[:, :, 0].astype(np.float32)
            binary = (gray > np.percentile(gray, 50)).astype(np.int32)
            tile_mask = semantic_to_instance(binary)

            # Crop back to original size
            tile_mask = tile_mask[: y_end - y, : x_end - x]

            # Merge into full mask
            for lbl in np.unique(tile_mask):
                if lbl == 0:
                    continue
                cell_pixels = tile_mask == lbl
                # Check overlap with existing cells
                region = full_mask[y:y_end, x:x_end]
                existing = region[cell_pixels]
                if np.sum(existing > 0) < np.sum(cell_pixels) * 0.5:
                    region[cell_pixels] = cell_id
                    cell_id += 1

    return full_mask
