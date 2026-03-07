"""Strategy 4: SegFormer-based segmentation with trainable N-channel input.

Uses a modified SegFormer (semantic segmentation) with variable input channels
to handle full IMC multi-channel stacks. Trained on datasets with ground truth
masks, then applied to all ROIs.

Post-processing converts semantic predictions to instance masks via distance
transform + watershed.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

__all__ = [
    "IMCSegmentationDataset",
    "build_segformer_model",
    "train_segformer",
    "predict_segformer",
    "postprocess_semantic_to_instance",
]

logger = logging.getLogger(__name__)


class IMCSegmentationDataset:
    """PyTorch Dataset for IMC segmentation training.

    Uses ``IMCPanelMapper`` for channel recognition and selects a fixed
    set of channels (DNA + key markers) across ROIs.

    Parameters
    ----------
    catalog
        DataFrame with tiff_path, channel_csv_path, mask_path columns.
    channels
        Channel names to select (resolved via IMCPanelMapper).
    crop_size
        Random crop size for training.
    augment
        Whether to apply data augmentation.
    """

    def __init__(
        self,
        catalog: pd.DataFrame,
        channels: list[str] | None = None,
        crop_size: int = 256,
        augment: bool = True,
    ):
        self.catalog = catalog[catalog["has_gt"]].reset_index(drop=True)
        self.channels = channels or ["DNA1", "DNA2", "CD3", "CD20", "PanCK", "SMA", "CD68"]
        self.crop_size = crop_size
        self.augment = augment

        if len(self.catalog) == 0:
            logger.warning("No ROIs with ground truth masks found in catalog")

    def __len__(self) -> int:
        return len(self.catalog)

    def __getitem__(self, idx: int) -> tuple[np.ndarray, np.ndarray]:
        """Return (image, mask) pair.

        Returns
        -------
        image : (n_channels, crop_size, crop_size) float32
        mask : (crop_size, crop_size) int64 — 0=bg, 1=cell, 2=boundary
        """
        row = self.catalog.iloc[idx]

        import tifffile

        from sc_tools.data.imc.benchmark.prepare import normalize_imc_intensity

        # Load full TIFF
        full_image = tifffile.imread(str(row["tiff_path"]))

        # Select channels via panel mapper
        selected = self._select_channels(full_image, row.get("channel_csv_path"))

        # Normalize: z-score per channel
        selected = normalize_imc_intensity(selected, method="zscore")

        # Load mask
        from sc_tools.bm.mask_io import load_mask

        mask = load_mask(row["mask_path"])

        # Create 3-class semantic mask: 0=bg, 1=cell interior, 2=boundary
        from skimage.segmentation import find_boundaries

        boundary = find_boundaries(mask, mode="inner")
        semantic = np.zeros_like(mask, dtype=np.int64)
        semantic[mask > 0] = 1
        semantic[boundary] = 2

        # Random crop
        selected, semantic = self._random_crop(selected, semantic)

        # Augmentation
        if self.augment:
            selected, semantic = self._augment(selected, semantic)

        return selected.astype(np.float32), semantic

    def _select_channels(
        self,
        full_image: np.ndarray,
        channel_csv_path: str | None,
    ) -> np.ndarray:
        """Select and reorder channels using IMCPanelMapper."""
        from sc_tools.ingest.imc import IMCPanelMapper

        mapper = IMCPanelMapper()
        if channel_csv_path and Path(channel_csv_path).is_file():
            mapper.from_full_csv(channel_csv_path)

        indices = []
        for ch_name in self.channels:
            idx = mapper.resolve(ch_name)
            if idx is not None and idx < full_image.shape[0]:
                indices.append(idx)
            else:
                # Use zeros for missing channels
                indices.append(-1)

        result = np.zeros(
            (len(self.channels), full_image.shape[1], full_image.shape[2]),
            dtype=np.float32,
        )
        for i, idx in enumerate(indices):
            if idx >= 0:
                result[i] = full_image[idx].astype(np.float32)

        return result

    def _random_crop(
        self,
        image: np.ndarray,
        mask: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Random crop of (C, H, W) image and (H, W) mask."""
        _, h, w = image.shape
        cs = self.crop_size

        if h <= cs or w <= cs:
            # Pad if needed
            pad_h = max(0, cs - h)
            pad_w = max(0, cs - w)
            image = np.pad(image, ((0, 0), (0, pad_h), (0, pad_w)), mode="reflect")
            mask = np.pad(mask, ((0, pad_h), (0, pad_w)), mode="reflect")
            h, w = image.shape[1], image.shape[2]

        y = np.random.randint(0, h - cs + 1)
        x = np.random.randint(0, w - cs + 1)

        return image[:, y : y + cs, x : x + cs], mask[y : y + cs, x : x + cs]

    def _augment(
        self,
        image: np.ndarray,
        mask: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Basic augmentation: random flip and rotation."""
        # Random horizontal flip
        if np.random.random() > 0.5:
            image = np.flip(image, axis=2).copy()
            mask = np.flip(mask, axis=1).copy()

        # Random vertical flip
        if np.random.random() > 0.5:
            image = np.flip(image, axis=1).copy()
            mask = np.flip(mask, axis=0).copy()

        # Random 90-degree rotation
        k = np.random.randint(0, 4)
        if k > 0:
            image = np.rot90(image, k, axes=(1, 2)).copy()
            mask = np.rot90(mask, k, axes=(0, 1)).copy()

        return image, mask


def build_segformer_model(
    n_input_channels: int = 7,
    n_classes: int = 3,
    model_size: str = "b0",
):
    """Build a modified SegFormer with N-channel input.

    Parameters
    ----------
    n_input_channels
        Number of input channels (e.g., 7 for DNA1+DNA2+CD3+CD20+PanCK+SMA+CD68).
    n_classes
        Number of output classes (3: background, cell, boundary).
    model_size
        SegFormer variant: ``"b0"`` through ``"b5"``.

    Returns
    -------
    PyTorch model.
    """
    try:
        import segmentation_models_pytorch as smp
    except ImportError as e:
        raise ImportError(
            "segmentation-models-pytorch required for SegFormer. "
            "Install with: pip install 'sc-tools[benchmark-extended]'"
        ) from e

    model = smp.create_model(
        "segformer",
        encoder_name=f"mit_{model_size}",
        in_channels=n_input_channels,
        classes=n_classes,
        encoder_weights=None,  # Train from scratch for N-channel
    )

    return model


def train_segformer(
    catalog: pd.DataFrame,
    config,
    output_dir: str | Path = "segformer_model",
    channels: list[str] | None = None,
) -> Path:
    """Train a SegFormer model on ROIs with ground truth masks.

    Parameters
    ----------
    catalog
        ROI catalog (must have has_gt=True rows).
    config
        BenchmarkConfig with vit_* settings.
    output_dir
        Where to save model checkpoints.
    channels
        Channel names to use.

    Returns
    -------
    Path to the best model checkpoint.
    """
    try:
        import torch
        import torch.nn as nn
        from torch.optim import AdamW
        from torch.optim.lr_scheduler import CosineAnnealingLR
        from torch.utils.data import DataLoader
    except ImportError as e:
        raise ImportError(
            "torch required for SegFormer training. "
            "Install with: pip install 'sc-tools[benchmark-extended]'"
        ) from e

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Split train/val
    gt_catalog = catalog[catalog["has_gt"]].reset_index(drop=True)
    n_train = int(len(gt_catalog) * config.vit_train_split)
    train_catalog = gt_catalog.iloc[:n_train]
    val_catalog = gt_catalog.iloc[n_train:]

    logger.info("Training: %d ROIs, Validation: %d ROIs", len(train_catalog), len(val_catalog))

    # Build datasets
    default_channels = channels or ["DNA1", "DNA2", "CD3", "CD20", "PanCK", "SMA", "CD68"]
    train_ds = IMCSegmentationDataset(train_catalog, channels=default_channels, augment=True)
    val_ds = IMCSegmentationDataset(val_catalog, channels=default_channels, augment=False)

    train_loader = DataLoader(train_ds, batch_size=config.batch_size, shuffle=True, num_workers=0)
    val_loader = DataLoader(val_ds, batch_size=config.batch_size, shuffle=False, num_workers=0)

    # Build model
    n_channels = len(default_channels)
    model = build_segformer_model(
        n_input_channels=n_channels, n_classes=3, model_size=config.vit_model_size
    )

    device = "cuda" if config.gpu and torch.cuda.is_available() else "cpu"
    model = model.to(device)

    # Loss: CE + Dice
    ce_loss = nn.CrossEntropyLoss()
    optimizer = AdamW(model.parameters(), lr=1e-4, weight_decay=0.01)
    scheduler = CosineAnnealingLR(optimizer, T_max=config.vit_n_epochs)

    # Training loop
    best_val_loss = float("inf")
    best_path = output_dir / "best_model.pt"

    for epoch in range(config.vit_n_epochs):
        # Train
        model.train()
        train_loss = 0.0
        n_batches = 0
        for images, masks in train_loader:
            images = torch.tensor(images, dtype=torch.float32).to(device)
            masks = torch.tensor(masks, dtype=torch.long).to(device)

            optimizer.zero_grad()
            outputs = model(images)

            # Handle different output formats
            if isinstance(outputs, dict):
                logits = outputs.get("logits", outputs.get("out", None))
                if logits is None:
                    logits = list(outputs.values())[0]
            else:
                logits = outputs

            # Resize logits to mask size if needed
            if logits.shape[-2:] != masks.shape[-2:]:
                logits = torch.nn.functional.interpolate(
                    logits, size=masks.shape[-2:], mode="bilinear", align_corners=False
                )

            loss = ce_loss(logits, masks) + _dice_loss(logits, masks, n_classes=3)
            loss.backward()
            optimizer.step()

            train_loss += loss.item()
            n_batches += 1

        scheduler.step()
        avg_train = train_loss / max(n_batches, 1)

        # Validate
        model.eval()
        val_loss = 0.0
        n_val = 0
        with torch.no_grad():
            for images, masks in val_loader:
                images = torch.tensor(images, dtype=torch.float32).to(device)
                masks = torch.tensor(masks, dtype=torch.long).to(device)

                outputs = model(images)
                logits = outputs if not isinstance(outputs, dict) else list(outputs.values())[0]

                if logits.shape[-2:] != masks.shape[-2:]:
                    logits = torch.nn.functional.interpolate(
                        logits, size=masks.shape[-2:], mode="bilinear", align_corners=False
                    )

                loss = ce_loss(logits, masks)
                val_loss += loss.item()
                n_val += 1

        avg_val = val_loss / max(n_val, 1)

        if avg_val < best_val_loss:
            best_val_loss = avg_val
            torch.save(model.state_dict(), best_path)

        if (epoch + 1) % 10 == 0 or epoch == 0:
            logger.info(
                "Epoch %d/%d: train_loss=%.4f, val_loss=%.4f (best=%.4f)",
                epoch + 1,
                config.vit_n_epochs,
                avg_train,
                avg_val,
                best_val_loss,
            )

    logger.info("Training complete. Best model saved to %s", best_path)
    return best_path


def predict_segformer(
    model_path: str | Path,
    tiff_path: str | Path,
    channel_csv_path: str | Path | None = None,
    channels: list[str] | None = None,
    tile_size: int = 256,
    overlap: int = 32,
    gpu: bool = False,
) -> np.ndarray:
    """Run SegFormer inference with sliding window.

    Parameters
    ----------
    model_path
        Path to trained model checkpoint.
    tiff_path
        Path to multi-channel TIFF.
    channel_csv_path
        Path to channel CSV.
    channels
        Channel names (must match training).
    tile_size
        Tile size for sliding window.
    overlap
        Overlap between tiles.
    gpu
        Whether to use GPU.

    Returns
    -------
    Instance segmentation mask (H, W), dtype uint32.
    """
    try:
        import torch
    except ImportError as e:
        raise ImportError("torch required: pip install 'sc-tools[benchmark-extended]'") from e

    import tifffile

    from sc_tools.data.imc.benchmark.prepare import normalize_imc_intensity

    default_channels = channels or ["DNA1", "DNA2", "CD3", "CD20", "PanCK", "SMA", "CD68"]

    # Load and select channels
    full_image = tifffile.imread(str(tiff_path))
    ds = IMCSegmentationDataset.__new__(IMCSegmentationDataset)
    ds.channels = default_channels
    selected = ds._select_channels(full_image, channel_csv_path)
    selected = normalize_imc_intensity(selected, method="zscore")

    # Load model
    device = "cuda" if gpu and torch.cuda.is_available() else "cpu"
    model = build_segformer_model(n_input_channels=len(default_channels), n_classes=3)
    model.load_state_dict(torch.load(str(model_path), map_location=device, weights_only=True))
    model = model.to(device)
    model.eval()

    # Sliding window inference
    _, h, w = selected.shape
    semantic = np.zeros((h, w), dtype=np.float32)
    counts = np.zeros((h, w), dtype=np.float32)

    step = tile_size - overlap
    with torch.no_grad():
        for y in range(0, h, step):
            for x in range(0, w, step):
                y_end = min(y + tile_size, h)
                x_end = min(x + tile_size, w)

                tile = selected[:, y:y_end, x:x_end]

                # Pad if needed
                pad_h = tile_size - tile.shape[1]
                pad_w = tile_size - tile.shape[2]
                if pad_h > 0 or pad_w > 0:
                    tile = np.pad(tile, ((0, 0), (0, pad_h), (0, pad_w)), mode="reflect")

                tensor = torch.tensor(tile[np.newaxis], dtype=torch.float32).to(device)
                output = model(tensor)

                if isinstance(output, dict):
                    logits = list(output.values())[0]
                else:
                    logits = output

                if logits.shape[-2:] != (tile_size, tile_size):
                    logits = torch.nn.functional.interpolate(
                        logits, size=(tile_size, tile_size), mode="bilinear", align_corners=False
                    )

                pred = logits[0, 1].cpu().numpy()  # cell class probability
                pred = pred[: y_end - y, : x_end - x]

                semantic[y:y_end, x:x_end] += pred
                counts[y:y_end, x:x_end] += 1.0

    # Average overlapping predictions
    semantic = semantic / np.maximum(counts, 1)

    # Post-process to instance mask
    return postprocess_semantic_to_instance(semantic)


def postprocess_semantic_to_instance(
    semantic_prob: np.ndarray,
    threshold: float = 0.5,
    min_area: int = 10,
) -> np.ndarray:
    """Convert semantic cell probability to instance mask.

    Parameters
    ----------
    semantic_prob
        2D probability map for cell class.
    threshold
        Threshold for cell/background.
    min_area
        Minimum cell area.

    Returns
    -------
    Instance mask (H, W), dtype uint32.
    """
    from sc_tools.bm.postprocess import filter_masks_by_area, semantic_to_instance

    binary = (semantic_prob > threshold).astype(np.int32)
    instances = semantic_to_instance(binary)
    instances = filter_masks_by_area(instances, min_area=min_area)

    return instances.astype(np.uint32)


def _dice_loss(logits, targets, n_classes: int = 3, smooth: float = 1.0):
    """Compute soft Dice loss for multi-class segmentation."""
    import torch.nn.functional as F

    probs = F.softmax(logits, dim=1)
    targets_one_hot = F.one_hot(targets, n_classes).permute(0, 3, 1, 2).float()

    dims = (0, 2, 3)
    intersection = (probs * targets_one_hot).sum(dim=dims)
    cardinality = (probs + targets_one_hot).sum(dim=dims)

    dice = (2.0 * intersection + smooth) / (cardinality + smooth)
    return 1.0 - dice.mean()
