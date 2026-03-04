"""Tests for sc_tools.bm.segmentation — cell segmentation quality metrics."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from sc_tools.bm.segmentation import (
    compare_segmentations,
    compute_detection_metrics,
    compute_marker_quality,
    compute_morphology_metrics,
    compute_segmentation_accuracy,
    compute_size_distribution,
    compute_spatial_coherence,
    score_segmentation,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _make_circular_mask(n_cells: int = 10, size: int = 200, radius: int = 8) -> np.ndarray:
    """Create a mask with roughly circular cells (known circularity ~1.0)."""
    mask = np.zeros((size, size), dtype=np.int32)
    rng = np.random.RandomState(42)
    cell_id = 1
    for _ in range(n_cells):
        cy = rng.randint(radius + 2, size - radius - 2)
        cx = rng.randint(radius + 2, size - radius - 2)
        yy, xx = np.ogrid[-cy : size - cy, -cx : size - cx]
        circle = (yy**2 + xx**2) <= radius**2
        if np.any(mask[circle] > 0):
            continue
        mask[circle] = cell_id
        cell_id += 1
    return mask


def _make_noisy_mask(base: np.ndarray, noise_frac: float = 0.1) -> np.ndarray:
    """Add random noise to a mask (merge/split cells)."""
    noisy = base.copy()
    rng = np.random.RandomState(99)
    n_pixels = int(noise_frac * np.sum(base > 0))
    ys, xs = np.where(base > 0)
    indices = rng.choice(len(ys), size=min(n_pixels, len(ys)), replace=False)
    for idx in indices:
        noisy[ys[idx], xs[idx]] = rng.randint(1, base.max() + 2)
    return noisy


def _make_intensity_image(
    mask: np.ndarray, n_markers: int = 3, signal_mean: float = 100.0, bg_mean: float = 10.0
) -> np.ndarray:
    """Create an intensity image with bright cells and dim background."""
    rng = np.random.RandomState(42)
    h, w = mask.shape
    img = np.zeros((h, w, n_markers), dtype=np.float32)
    for ch in range(n_markers):
        img[:, :, ch] = rng.poisson(bg_mean, size=(h, w)).astype(np.float32)
        img[:, :, ch][mask > 0] += rng.poisson(signal_mean, size=np.sum(mask > 0)).astype(
            np.float32
        )
    return img


# ---------------------------------------------------------------------------
# Morphology metrics
# ---------------------------------------------------------------------------


class TestMorphologyMetrics:
    def test_returns_dataframe(self):
        mask = _make_circular_mask(5)
        df = compute_morphology_metrics(mask)
        assert isinstance(df, pd.DataFrame)
        assert set(df.columns) >= {
            "label",
            "area",
            "perimeter",
            "circularity",
            "solidity",
            "eccentricity",
            "aspect_ratio",
        }

    def test_correct_cell_count(self):
        mask = _make_circular_mask(8)
        df = compute_morphology_metrics(mask)
        n_cells = len(np.unique(mask)) - 1  # exclude bg
        assert len(df) == n_cells

    def test_circular_cells_high_circularity(self):
        mask = _make_circular_mask(10, radius=12)
        df = compute_morphology_metrics(mask)
        # Circular cells should have circularity close to 1.0
        median_circ = df["circularity"].median()
        assert median_circ > 0.8, f"Expected high circularity, got {median_circ}"

    def test_empty_mask(self):
        mask = np.zeros((50, 50), dtype=np.int32)
        df = compute_morphology_metrics(mask)
        assert len(df) == 0


# ---------------------------------------------------------------------------
# Marker quality
# ---------------------------------------------------------------------------


class TestMarkerQuality:
    def test_returns_dataframe(self):
        mask = _make_circular_mask(15)
        img = _make_intensity_image(mask, n_markers=2)
        df = compute_marker_quality(mask, img, marker_names=["CD3", "CD20"])
        assert isinstance(df, pd.DataFrame)
        assert set(df.columns) >= {"marker", "snr", "mean_positive", "mean_negative", "purity"}
        assert len(df) == 2

    def test_positive_snr(self):
        mask = _make_circular_mask(20)
        img = _make_intensity_image(mask, signal_mean=200, bg_mean=5)
        df = compute_marker_quality(mask, img)
        # Strong signal should give positive SNR
        assert all(df["snr"] > 0), f"Expected positive SNR, got {df['snr'].tolist()}"

    def test_single_channel(self):
        mask = _make_circular_mask(10)
        img = _make_intensity_image(mask, n_markers=1)[:, :, 0]
        df = compute_marker_quality(mask, img)
        assert len(df) == 1


# ---------------------------------------------------------------------------
# Spatial coherence
# ---------------------------------------------------------------------------


class TestSpatialCoherence:
    def test_returns_dict(self):
        mask = _make_circular_mask(20)
        result = compute_spatial_coherence(mask)
        assert isinstance(result, dict)
        assert set(result.keys()) >= {
            "cell_density_cv",
            "neighbor_count_mean",
            "neighbor_count_std",
            "boundary_regularity",
        }

    def test_few_cells_handled(self):
        mask = np.zeros((50, 50), dtype=np.int32)
        mask[10:15, 10:15] = 1
        result = compute_spatial_coherence(mask)
        assert result["neighbor_count_mean"] == 0.0


# ---------------------------------------------------------------------------
# Size distribution
# ---------------------------------------------------------------------------


class TestSizeDistribution:
    def test_returns_dict(self):
        mask = _make_circular_mask(10)
        result = compute_size_distribution(mask)
        assert isinstance(result, dict)
        assert "area_cv" in result
        assert "pct_outlier" in result
        assert "median_area" in result
        assert "n_cells" in result

    def test_cell_count_matches(self):
        mask = _make_circular_mask(8)
        result = compute_size_distribution(mask)
        n_cells = len(np.unique(mask)) - 1
        assert result["n_cells"] == n_cells

    def test_uniform_cells_low_outlier(self):
        mask = _make_circular_mask(15, radius=10)
        result = compute_size_distribution(mask)
        assert result["pct_outlier"] < 30, f"Expected low outlier %, got {result['pct_outlier']}"

    def test_empty_mask(self):
        result = compute_size_distribution(np.zeros((50, 50), dtype=np.int32))
        assert result["n_cells"] == 0


# ---------------------------------------------------------------------------
# Detection metrics (with ground truth)
# ---------------------------------------------------------------------------


class TestDetectionMetrics:
    def test_perfect_detection(self):
        mask = _make_circular_mask(5)
        result = compute_detection_metrics(mask, mask)
        assert result["precision"] == 1.0
        assert result["recall"] == 1.0
        assert result["f1"] == 1.0

    def test_no_predictions(self):
        gt = _make_circular_mask(5)
        pred = np.zeros_like(gt)
        result = compute_detection_metrics(pred, gt)
        assert result["f1"] == 0.0
        assert result["n_fn"] > 0

    def test_no_gt(self):
        pred = _make_circular_mask(5)
        gt = np.zeros_like(pred)
        result = compute_detection_metrics(pred, gt)
        assert result["f1"] == 0.0
        assert result["n_fp"] > 0

    def test_partial_detection(self):
        gt = _make_circular_mask(10)
        # Remove half the cells from pred
        pred = gt.copy()
        labels = np.unique(pred)
        labels = labels[labels > 0]
        for lbl in labels[: len(labels) // 2]:
            pred[pred == lbl] = 0
        result = compute_detection_metrics(pred, gt)
        assert 0 < result["f1"] < 1.0


# ---------------------------------------------------------------------------
# Segmentation accuracy
# ---------------------------------------------------------------------------


class TestSegmentationAccuracy:
    def test_perfect_accuracy(self):
        mask = _make_circular_mask(5)
        result = compute_segmentation_accuracy(mask, mask)
        assert result["mean_iou"] == pytest.approx(1.0, abs=0.01)
        assert result["mean_dice"] == pytest.approx(1.0, abs=0.01)
        assert result["ap_50"] == pytest.approx(1.0, abs=0.01)

    def test_empty_gives_zeros(self):
        result = compute_segmentation_accuracy(
            np.zeros((50, 50), dtype=np.int32),
            _make_circular_mask(5),
        )
        assert result["mean_iou"] == 0.0


# ---------------------------------------------------------------------------
# score_segmentation (all-in-one)
# ---------------------------------------------------------------------------


class TestScoreSegmentation:
    def test_basic_keys(self):
        mask = _make_circular_mask(10)
        result = score_segmentation(mask)
        assert "morphology" in result
        assert "spatial_coherence" in result
        assert "size_distribution" in result

    def test_with_intensity(self):
        mask = _make_circular_mask(10)
        img = _make_intensity_image(mask)
        result = score_segmentation(mask, intensity_image=img)
        assert "marker_quality" in result

    def test_with_gt(self):
        mask = _make_circular_mask(10)
        result = score_segmentation(mask, gt_mask=mask)
        assert "detection" in result
        assert "accuracy" in result


# ---------------------------------------------------------------------------
# compare_segmentations
# ---------------------------------------------------------------------------


class TestCompareSegmentations:
    def test_returns_sorted_df(self):
        clean = _make_circular_mask(15)
        noisy = _make_noisy_mask(clean, noise_frac=0.2)
        df = compare_segmentations({"clean": clean, "noisy": noisy})
        assert isinstance(df, pd.DataFrame)
        assert "composite_score" in df.columns
        assert len(df) == 2
        # Should be sorted descending
        assert df["composite_score"].iloc[0] >= df["composite_score"].iloc[1]

    def test_composite_score_range(self):
        masks = {
            "a": _make_circular_mask(10),
            "b": _make_circular_mask(10, radius=5),
            "c": _make_noisy_mask(_make_circular_mask(10)),
        }
        df = compare_segmentations(masks)
        assert all(0 <= s <= 100 for s in df["composite_score"])

    def test_single_method(self):
        mask = _make_circular_mask(10)
        df = compare_segmentations({"only": mask})
        assert len(df) == 1
        assert df["composite_score"].iloc[0] == 50.0
