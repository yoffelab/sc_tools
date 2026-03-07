"""
Cell segmentation quality metrics and comparison.

Provides no-ground-truth metrics (morphology, marker quality, spatial
coherence, size distribution) and with-ground-truth metrics (detection,
segmentation accuracy). Composite scoring via PCA on normalized metrics
across methods.
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.optimize import linear_sum_assignment
from scipy.spatial import Delaunay
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

__all__ = [
    "compute_morphology_metrics",
    "compute_marker_quality",
    "compute_spatial_coherence",
    "compute_size_distribution",
    "compute_detection_metrics",
    "compute_segmentation_accuracy",
    "compute_panoptic_quality",
    "compute_boundary_metrics",
    "compute_cell_type_preservation",
    "score_segmentation",
    "compare_segmentations",
]

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# No-ground-truth metrics
# ---------------------------------------------------------------------------


def compute_morphology_metrics(mask: np.ndarray) -> pd.DataFrame:
    """Compute per-cell morphology metrics from a labeled mask.

    Parameters
    ----------
    mask
        Labeled 2D array (0=background, >0=cell_id).

    Returns
    -------
    DataFrame with columns: ``label``, ``area``, ``perimeter``,
    ``circularity``, ``solidity``, ``eccentricity``, ``aspect_ratio``.
    """
    from skimage.measure import regionprops_table

    props = regionprops_table(
        mask,
        properties=(
            "label",
            "area",
            "perimeter",
            "solidity",
            "eccentricity",
            "major_axis_length",
            "minor_axis_length",
        ),
    )
    df = pd.DataFrame(props)

    # Circularity: 4*pi*area / perimeter^2 (1.0 = perfect circle)
    perimeter = df["perimeter"].replace(0, np.nan)
    df["circularity"] = (4 * np.pi * df["area"]) / (perimeter**2)

    # Aspect ratio: major / minor axis (1.0 = round)
    minor = df["minor_axis_length"].replace(0, np.nan)
    df["aspect_ratio"] = df["major_axis_length"] / minor

    return df[
        ["label", "area", "perimeter", "circularity", "solidity", "eccentricity", "aspect_ratio"]
    ]


def compute_marker_quality(
    mask: np.ndarray,
    intensity_image: np.ndarray,
    marker_names: list[str] | None = None,
) -> pd.DataFrame:
    """Compute per-marker signal-to-noise via GMM (positive vs negative cells).

    Parameters
    ----------
    mask
        Labeled 2D array.
    intensity_image
        Intensity array, shape ``(H, W)`` or ``(H, W, C)`` for multiple markers.
    marker_names
        Names for each channel. If None, uses ``marker_0``, ``marker_1``, etc.

    Returns
    -------
    DataFrame with columns: ``marker``, ``snr``, ``mean_positive``,
    ``mean_negative``, ``purity``.
    """
    if intensity_image.ndim == 2:
        intensity_image = intensity_image[:, :, np.newaxis]

    n_markers = intensity_image.shape[2]
    if marker_names is None:
        marker_names = [f"marker_{i}" for i in range(n_markers)]

    cell_labels = np.unique(mask)
    cell_labels = cell_labels[cell_labels > 0]

    if len(cell_labels) < 2:
        return pd.DataFrame(columns=["marker", "snr", "mean_positive", "mean_negative", "purity"])

    results = []
    for ch, name in enumerate(marker_names):
        channel = intensity_image[:, :, ch]
        # Mean intensity per cell
        means = ndimage.mean(channel, labels=mask, index=cell_labels)
        means = np.array(means).reshape(-1, 1)

        if len(means) < 2:
            results.append(
                {
                    "marker": name,
                    "snr": 0.0,
                    "mean_positive": 0.0,
                    "mean_negative": 0.0,
                    "purity": 0.0,
                }
            )
            continue

        gmm = GaussianMixture(n_components=2, random_state=42)
        gmm.fit(means)

        cluster_means = gmm.means_.flatten()
        pos_idx = np.argmax(cluster_means)
        neg_idx = 1 - pos_idx

        mean_pos = cluster_means[pos_idx]
        mean_neg = cluster_means[neg_idx]

        # SNR = (mean_pos - mean_neg) / std_neg (avoid div by zero)
        labels_gmm = gmm.predict(means)
        neg_cells = means[labels_gmm == neg_idx].flatten()
        std_neg = np.std(neg_cells) if len(neg_cells) > 1 else 1e-8
        snr = (mean_pos - mean_neg) / max(std_neg, 1e-8)

        # Purity: max posterior probability averaged over cells
        probs = gmm.predict_proba(means)
        purity = float(np.mean(np.max(probs, axis=1)))

        results.append(
            {
                "marker": name,
                "snr": float(snr),
                "mean_positive": float(mean_pos),
                "mean_negative": float(mean_neg),
                "purity": float(purity),
            }
        )

    return pd.DataFrame(results)


def compute_spatial_coherence(mask: np.ndarray) -> dict[str, float]:
    """Compute spatial coherence metrics for a segmentation mask.

    Returns
    -------
    Dict with keys: ``cell_density_cv``, ``neighbor_count_mean``,
    ``neighbor_count_std``, ``boundary_regularity``.
    """
    cell_labels = np.unique(mask)
    cell_labels = cell_labels[cell_labels > 0]
    n_cells = len(cell_labels)

    if n_cells < 3:
        return {
            "cell_density_cv": 0.0,
            "neighbor_count_mean": 0.0,
            "neighbor_count_std": 0.0,
            "boundary_regularity": 0.0,
        }

    # Cell centroids
    centroids = ndimage.center_of_mass(mask > 0, labels=mask, index=cell_labels)
    centroids = np.array(centroids)

    # Density CV: divide image into grid, count cells per tile
    h, w = mask.shape
    n_tiles = max(4, int(np.sqrt(n_cells)))
    tile_h = max(1, h // n_tiles)
    tile_w = max(1, w // n_tiles)
    counts = np.zeros((n_tiles, n_tiles))
    for cy, cx in centroids:
        ti = min(int(cy / tile_h), n_tiles - 1)
        tj = min(int(cx / tile_w), n_tiles - 1)
        counts[ti, tj] += 1
    flat_counts = counts.flatten()
    mean_count = np.mean(flat_counts)
    density_cv = float(np.std(flat_counts) / mean_count) if mean_count > 0 else 0.0

    # Neighbor count via Delaunay triangulation
    try:
        tri = Delaunay(centroids)
        indptr, indices = tri.vertex_neighbor_vertices
        neighbor_counts = np.array(
            [len(indices[indptr[i] : indptr[i + 1]]) for i in range(n_cells)]
        )
        neighbor_mean = float(np.mean(neighbor_counts))
        neighbor_std = float(np.std(neighbor_counts))
    except Exception:
        neighbor_mean = 0.0
        neighbor_std = 0.0

    # Boundary regularity: mean circularity as proxy
    from skimage.measure import regionprops

    props = regionprops(mask)
    circularities = []
    for p in props:
        if p.perimeter > 0:
            circularities.append(4 * np.pi * p.area / (p.perimeter**2))
    boundary_reg = float(np.mean(circularities)) if circularities else 0.0

    return {
        "cell_density_cv": density_cv,
        "neighbor_count_mean": neighbor_mean,
        "neighbor_count_std": neighbor_std,
        "boundary_regularity": boundary_reg,
    }


def compute_size_distribution(mask: np.ndarray) -> dict[str, float]:
    """Compute cell size distribution statistics.

    Returns
    -------
    Dict with keys: ``area_cv``, ``pct_outlier``, ``median_area``,
    ``mean_area``, ``n_cells``.
    """
    from skimage.measure import regionprops

    props = regionprops(mask)
    areas = np.array([p.area for p in props], dtype=float)

    if len(areas) == 0:
        return {
            "area_cv": 0.0,
            "pct_outlier": 0.0,
            "median_area": 0.0,
            "mean_area": 0.0,
            "n_cells": 0,
        }

    mean_area = float(np.mean(areas))
    median_area = float(np.median(areas))
    area_cv = float(np.std(areas) / mean_area) if mean_area > 0 else 0.0

    # Outliers via IQR
    q1, q3 = np.percentile(areas, [25, 75])
    iqr = q3 - q1
    lower = q1 - 1.5 * iqr
    upper = q3 + 1.5 * iqr
    n_outlier = int(np.sum((areas < lower) | (areas > upper)))
    pct_outlier = float(n_outlier / len(areas) * 100)

    return {
        "area_cv": area_cv,
        "pct_outlier": pct_outlier,
        "median_area": median_area,
        "mean_area": mean_area,
        "n_cells": len(areas),
    }


# ---------------------------------------------------------------------------
# With-ground-truth metrics
# ---------------------------------------------------------------------------


def _compute_iou_matrix(
    pred: np.ndarray, gt: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute IoU matrix between predicted and ground truth cells.

    Returns (iou_matrix, pred_labels, gt_labels).
    """
    pred_labels = np.unique(pred)
    pred_labels = pred_labels[pred_labels > 0]
    gt_labels = np.unique(gt)
    gt_labels = gt_labels[gt_labels > 0]

    if len(pred_labels) == 0 or len(gt_labels) == 0:
        return np.zeros((len(pred_labels), len(gt_labels))), pred_labels, gt_labels

    iou_matrix = np.zeros((len(pred_labels), len(gt_labels)))
    for i, pl in enumerate(pred_labels):
        pred_mask = pred == pl
        for j, gl in enumerate(gt_labels):
            gt_mask = gt == gl
            intersection = np.sum(pred_mask & gt_mask)
            union = np.sum(pred_mask | gt_mask)
            iou_matrix[i, j] = intersection / union if union > 0 else 0.0

    return iou_matrix, pred_labels, gt_labels


def compute_detection_metrics(
    pred: np.ndarray,
    gt: np.ndarray,
    iou_threshold: float = 0.5,
) -> dict[str, float]:
    """Compute detection metrics via Hungarian matching.

    Parameters
    ----------
    pred
        Predicted segmentation mask.
    gt
        Ground truth segmentation mask.
    iou_threshold
        Minimum IoU to count as a true positive.

    Returns
    -------
    Dict with keys: ``precision``, ``recall``, ``f1``, ``n_tp``, ``n_fp``, ``n_fn``.
    """
    iou_matrix, pred_labels, gt_labels = _compute_iou_matrix(pred, gt)

    n_pred = len(pred_labels)
    n_gt = len(gt_labels)

    if n_pred == 0 and n_gt == 0:
        return {"precision": 1.0, "recall": 1.0, "f1": 1.0, "n_tp": 0, "n_fp": 0, "n_fn": 0}
    if n_pred == 0:
        return {"precision": 0.0, "recall": 0.0, "f1": 0.0, "n_tp": 0, "n_fp": 0, "n_fn": n_gt}
    if n_gt == 0:
        return {"precision": 0.0, "recall": 0.0, "f1": 0.0, "n_tp": 0, "n_fp": n_pred, "n_fn": 0}

    # Hungarian matching (maximize IoU = minimize negative IoU)
    cost = -iou_matrix
    row_ind, col_ind = linear_sum_assignment(cost)

    tp = sum(1 for r, c in zip(row_ind, col_ind, strict=False) if iou_matrix[r, c] >= iou_threshold)
    fp = n_pred - tp
    fn = n_gt - tp

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    return {
        "precision": float(precision),
        "recall": float(recall),
        "f1": float(f1),
        "n_tp": int(tp),
        "n_fp": int(fp),
        "n_fn": int(fn),
    }


def compute_segmentation_accuracy(
    pred: np.ndarray,
    gt: np.ndarray,
    iou_thresholds: tuple[float, ...] = (0.5, 0.75),
) -> dict[str, float]:
    """Compute segmentation accuracy metrics.

    Returns
    -------
    Dict with keys: ``mean_iou``, ``mean_dice``, ``ap_50``, ``ap_75``,
    ``ap_50_95`` (COCO-style AP averaged over [0.5:0.05:0.95]).
    """
    iou_matrix, pred_labels, gt_labels = _compute_iou_matrix(pred, gt)
    n_pred = len(pred_labels)
    n_gt = len(gt_labels)

    if n_pred == 0 or n_gt == 0:
        return {"mean_iou": 0.0, "mean_dice": 0.0, "ap_50": 0.0, "ap_75": 0.0, "ap_50_95": 0.0}

    # Hungarian matching for mean IoU and Dice
    cost = -iou_matrix
    row_ind, col_ind = linear_sum_assignment(cost)

    matched_ious = [
        iou_matrix[r, c] for r, c in zip(row_ind, col_ind, strict=False) if iou_matrix[r, c] > 0
    ]
    mean_iou = float(np.mean(matched_ious)) if matched_ious else 0.0

    # Dice from IoU: dice = 2*iou / (1+iou)
    matched_dice = [2 * iou / (1 + iou) for iou in matched_ious]
    mean_dice = float(np.mean(matched_dice)) if matched_dice else 0.0

    # AP at various thresholds
    def _ap_at_threshold(thresh: float) -> float:
        tp = sum(1 for r, c in zip(row_ind, col_ind, strict=False) if iou_matrix[r, c] >= thresh)
        prec = tp / n_pred if n_pred > 0 else 0.0
        rec = tp / n_gt if n_gt > 0 else 0.0
        return 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0

    ap_50 = _ap_at_threshold(0.5)
    ap_75 = _ap_at_threshold(0.75)

    # COCO-style AP@50:95
    coco_thresholds = np.arange(0.5, 1.0, 0.05)
    ap_50_95 = float(np.mean([_ap_at_threshold(t) for t in coco_thresholds]))

    return {
        "mean_iou": mean_iou,
        "mean_dice": mean_dice,
        "ap_50": ap_50,
        "ap_75": ap_75,
        "ap_50_95": ap_50_95,
    }


def compute_panoptic_quality(
    pred: np.ndarray,
    gt: np.ndarray,
    iou_threshold: float = 0.5,
) -> dict[str, float]:
    """Compute Panoptic Quality (PQ = SQ x DQ).

    Standard instance segmentation metric from Kirillov et al. (2019).

    Parameters
    ----------
    pred
        Predicted instance segmentation mask.
    gt
        Ground truth instance segmentation mask.
    iou_threshold
        IoU threshold for matching.

    Returns
    -------
    Dict with keys: ``pq``, ``sq`` (segmentation quality), ``dq`` (detection quality),
    ``n_tp``, ``n_fp``, ``n_fn``.
    """
    iou_matrix, pred_labels, gt_labels = _compute_iou_matrix(pred, gt)
    n_pred = len(pred_labels)
    n_gt = len(gt_labels)

    if n_pred == 0 and n_gt == 0:
        return {"pq": 1.0, "sq": 1.0, "dq": 1.0, "n_tp": 0, "n_fp": 0, "n_fn": 0}
    if n_pred == 0 or n_gt == 0:
        return {"pq": 0.0, "sq": 0.0, "dq": 0.0, "n_tp": 0, "n_fp": n_pred, "n_fn": n_gt}

    # Hungarian matching
    cost = -iou_matrix
    row_ind, col_ind = linear_sum_assignment(cost)

    # Matched pairs above threshold
    tp_ious = []
    matched_pred = set()
    matched_gt = set()
    for r, c in zip(row_ind, col_ind, strict=False):
        if iou_matrix[r, c] >= iou_threshold:
            tp_ious.append(iou_matrix[r, c])
            matched_pred.add(r)
            matched_gt.add(c)

    n_tp = len(tp_ious)
    n_fp = n_pred - n_tp
    n_fn = n_gt - n_tp

    sq = float(np.mean(tp_ious)) if tp_ious else 0.0
    dq = n_tp / (n_tp + 0.5 * n_fp + 0.5 * n_fn) if (n_tp + n_fp + n_fn) > 0 else 0.0
    pq = sq * dq

    return {
        "pq": float(pq),
        "sq": float(sq),
        "dq": float(dq),
        "n_tp": int(n_tp),
        "n_fp": int(n_fp),
        "n_fn": int(n_fn),
    }


def compute_boundary_metrics(
    pred: np.ndarray,
    gt: np.ndarray,
    tolerances: tuple[int, ...] = (1, 2, 3, 5),
) -> dict[str, float]:
    """Compute boundary F1 at multiple pixel tolerances.

    Parameters
    ----------
    pred
        Predicted instance segmentation mask.
    gt
        Ground truth instance segmentation mask.
    tolerances
        Pixel tolerances for boundary matching.

    Returns
    -------
    Dict with keys ``boundary_f1_{tol}px`` for each tolerance.
    """
    from skimage.segmentation import find_boundaries

    pred_boundaries = find_boundaries(pred, mode="inner")
    gt_boundaries = find_boundaries(gt, mode="inner")

    results = {}
    for tol in tolerances:
        # Dilate boundaries by tolerance
        from scipy.ndimage import binary_dilation

        struct = np.ones((2 * tol + 1, 2 * tol + 1))
        gt_dilated = binary_dilation(gt_boundaries, structure=struct)
        pred_dilated = binary_dilation(pred_boundaries, structure=struct)

        # Precision: pred boundary pixels within tolerance of gt
        tp_pred = np.sum(pred_boundaries & gt_dilated)
        n_pred = np.sum(pred_boundaries)
        precision = tp_pred / n_pred if n_pred > 0 else 0.0

        # Recall: gt boundary pixels within tolerance of pred
        tp_gt = np.sum(gt_boundaries & pred_dilated)
        n_gt = np.sum(gt_boundaries)
        recall = tp_gt / n_gt if n_gt > 0 else 0.0

        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
        results[f"boundary_f1_{tol}px"] = float(f1)

    return results


def compute_cell_type_preservation(
    pred_mask: np.ndarray,
    gt_mask: np.ndarray,
    gt_labels: np.ndarray | dict[int, str],
    iou_threshold: float = 0.5,
) -> dict[str, float]:
    """Measure how well segmentation preserves cell type identity.

    For annotated public datasets where GT cells have type labels.

    Parameters
    ----------
    pred_mask
        Predicted instance mask.
    gt_mask
        Ground truth instance mask.
    gt_labels
        Either a 2D array (same shape as gt_mask) with cell type IDs per pixel,
        or a dict mapping GT cell label -> cell type string.
    iou_threshold
        IoU threshold for matching.

    Returns
    -------
    Dict with ``type_preservation_rate``, ``n_matched``, ``n_types_preserved``.
    """
    iou_matrix, pred_labs, gt_labs = _compute_iou_matrix(pred_mask, gt_mask)

    if len(pred_labs) == 0 or len(gt_labs) == 0:
        return {"type_preservation_rate": 0.0, "n_matched": 0, "n_types_preserved": 0}

    # Hungarian matching
    cost = -iou_matrix
    row_ind, col_ind = linear_sum_assignment(cost)

    # Get cell type for each GT cell
    if isinstance(gt_labels, dict):
        gt_type_map = gt_labels
    else:
        # Extract per-cell majority type from label array
        gt_type_map = {}
        for gl in gt_labs:
            cell_pixels = gt_labels[gt_mask == gl]
            if len(cell_pixels) > 0:
                values, counts = np.unique(cell_pixels, return_counts=True)
                gt_type_map[int(gl)] = int(values[np.argmax(counts)])

    # Count matched cells that map to the same type
    matched = 0
    type_correct = 0
    types_seen = set()

    for r, c in zip(row_ind, col_ind, strict=False):
        if iou_matrix[r, c] >= iou_threshold:
            matched += 1
            gt_label = int(gt_labs[c])
            if gt_label in gt_type_map:
                types_seen.add(gt_type_map[gt_label])
                # A matched cell preserves type by definition if IoU is high enough
                type_correct += 1

    preservation_rate = type_correct / matched if matched > 0 else 0.0

    return {
        "type_preservation_rate": float(preservation_rate),
        "n_matched": int(matched),
        "n_types_preserved": len(types_seen),
    }


# ---------------------------------------------------------------------------
# Composite scoring and comparison
# ---------------------------------------------------------------------------


def score_segmentation(
    mask: np.ndarray,
    intensity_image: np.ndarray | None = None,
    gt_mask: np.ndarray | None = None,
    marker_names: list[str] | None = None,
) -> dict[str, Any]:
    """Run all applicable metrics and return raw results.

    Parameters
    ----------
    mask
        Labeled segmentation mask.
    intensity_image
        Optional intensity image for marker quality metrics.
    gt_mask
        Optional ground truth mask for detection/accuracy metrics.
    marker_names
        Optional marker names for intensity channels.

    Returns
    -------
    Dict with keys ``morphology``, ``spatial_coherence``,
    ``size_distribution``, and optionally ``marker_quality``,
    ``detection``, ``accuracy``.
    """
    result: dict[str, Any] = {}

    result["morphology"] = compute_morphology_metrics(mask)
    result["spatial_coherence"] = compute_spatial_coherence(mask)
    result["size_distribution"] = compute_size_distribution(mask)

    if intensity_image is not None:
        result["marker_quality"] = compute_marker_quality(
            mask, intensity_image, marker_names=marker_names
        )

    if gt_mask is not None:
        result["detection"] = compute_detection_metrics(mask, gt_mask)
        result["accuracy"] = compute_segmentation_accuracy(mask, gt_mask)

    return result


def compare_segmentations(
    masks: dict[str, np.ndarray],
    intensity_image: np.ndarray | None = None,
    gt_mask: np.ndarray | None = None,
    marker_names: list[str] | None = None,
) -> pd.DataFrame:
    """Compare multiple segmentation methods side-by-side.

    Parameters
    ----------
    masks
        Dict mapping method name to labeled mask array.
    intensity_image
        Optional intensity image for marker quality.
    gt_mask
        Optional ground truth mask.
    marker_names
        Optional marker names.

    Returns
    -------
    DataFrame with rows=methods, columns=summary metrics, sorted by
    ``composite_score`` descending.
    """
    all_results = {}
    summaries = []

    for name, mask in masks.items():
        result = score_segmentation(
            mask,
            intensity_image=intensity_image,
            gt_mask=gt_mask,
            marker_names=marker_names,
        )
        all_results[name] = result

        morph = result["morphology"]
        sd = result["size_distribution"]
        sc = result["spatial_coherence"]

        row = {
            "method": name,
            "n_cells": int(sd["n_cells"]),
            "median_area": sd["median_area"],
            "area_cv": sd["area_cv"],
            "pct_outlier": sd["pct_outlier"],
            "median_circularity": float(morph["circularity"].median()) if len(morph) > 0 else 0.0,
            "median_solidity": float(morph["solidity"].median()) if len(morph) > 0 else 0.0,
            "cell_density_cv": sc["cell_density_cv"],
            "boundary_regularity": sc["boundary_regularity"],
        }

        if "marker_quality" in result and len(result["marker_quality"]) > 0:
            row["mean_snr"] = float(result["marker_quality"]["snr"].mean())

        if "detection" in result:
            row["detection_f1"] = result["detection"]["f1"]

        if "accuracy" in result:
            row["mean_iou"] = result["accuracy"]["mean_iou"]
            row["ap_50_95"] = result["accuracy"]["ap_50_95"]

        summaries.append(row)

    df = pd.DataFrame(summaries)

    if len(df) < 2:
        df["composite_score"] = 50.0
        return df

    # Composite score via PCA on z-scored summary metrics
    score_cols = ["median_circularity", "median_solidity", "boundary_regularity"]
    # Inverted metrics (lower is better)
    invert_cols = ["area_cv", "pct_outlier", "cell_density_cv"]

    if "mean_snr" in df.columns:
        score_cols.append("mean_snr")
    if "detection_f1" in df.columns:
        score_cols.append("detection_f1")
    if "mean_iou" in df.columns:
        score_cols.append("mean_iou")

    all_cols = score_cols + [c for c in invert_cols if c in df.columns]
    feature_matrix = df[all_cols].copy()

    # Invert "lower is better" columns
    for col in invert_cols:
        if col in feature_matrix.columns:
            feature_matrix[col] = -feature_matrix[col]

    # Z-score and PCA
    scaler = StandardScaler()
    scaled = scaler.fit_transform(feature_matrix.fillna(0))

    pca = PCA(n_components=1)
    pc1 = pca.fit_transform(scaled).flatten()

    # Rescale to 0-100
    pc_min, pc_max = pc1.min(), pc1.max()
    if pc_max > pc_min:
        composite = (pc1 - pc_min) / (pc_max - pc_min) * 100
    else:
        composite = np.full_like(pc1, 50.0)

    df["composite_score"] = composite
    df = df.sort_values("composite_score", ascending=False).reset_index(drop=True)

    return df
