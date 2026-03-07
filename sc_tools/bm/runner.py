"""Benchmark orchestration: run segmentation strategies across ROIs.

Supports resume (skip completed ROIs), intermediate result saving,
and aggregation across datasets and tissue types.
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd

__all__ = [
    "run_single_roi",
    "run_benchmark",
    "aggregate_results",
]

logger = logging.getLogger(__name__)


def run_single_roi(
    tiff_path: str | Path,
    strategy: int,
    methods: list[str],
    channel_csv_path: str | Path | None = None,
    panel_mapper=None,
    prob_path: str | Path | None = None,
    gt_mask_path: str | Path | None = None,
    intensity_image: np.ndarray | None = None,
    gpu: bool = False,
    prob_map_method: str = "gaussian",
    output_dir: str | Path | None = None,
) -> list[dict]:
    """Run segmentation methods on a single ROI and evaluate.

    Parameters
    ----------
    tiff_path
        Path to multi-channel TIFF.
    strategy
        Strategy number (1-4).
    methods
        Method names to run.
    channel_csv_path
        Path to channel CSV.
    panel_mapper
        Optional IMCPanelMapper.
    prob_path
        Path to Ilastik probability map (Strategy 1 only).
    gt_mask_path
        Path to ground truth mask for evaluation.
    intensity_image
        Optional intensity image for marker quality metrics.
    gpu
        Whether to use GPU.
    prob_map_method
        Prob map generation method (Strategy 2).
    output_dir
        Where to save masks. None = do not save.

    Returns
    -------
    List of result dicts, one per method.
    """
    import tifffile

    from sc_tools.bm.segmentation import score_segmentation

    results = []

    # Load ground truth if available
    gt_mask = None
    if gt_mask_path is not None:
        gt_mask_path = Path(gt_mask_path)
        if gt_mask_path.is_file():
            from sc_tools.bm.mask_io import load_mask

            gt_mask = load_mask(gt_mask_path)

    # Load intensity image for marker quality
    if intensity_image is None:
        try:
            intensity_image = tifffile.imread(str(tiff_path))
        except Exception:
            pass

    # Get masks from the appropriate strategy
    masks = _run_strategy(
        strategy=strategy,
        tiff_path=tiff_path,
        methods=methods,
        channel_csv_path=channel_csv_path,
        panel_mapper=panel_mapper,
        prob_path=prob_path,
        gpu=gpu,
        prob_map_method=prob_map_method,
    )

    # Evaluate each mask
    for method_name, mask in masks.items():
        t0 = time.time()

        # Prepare intensity for marker quality (H, W, C format)
        int_img = None
        if intensity_image is not None and intensity_image.ndim == 3:
            # (C, H, W) -> (H, W, C) for marker quality
            int_img = np.moveaxis(intensity_image, 0, -1)

        scores = score_segmentation(
            mask,
            intensity_image=int_img,
            gt_mask=gt_mask,
        )

        elapsed = time.time() - t0

        row = {
            "strategy": strategy,
            "method": method_name,
            "n_cells": int(scores["size_distribution"]["n_cells"]),
            "median_area": scores["size_distribution"]["median_area"],
            "area_cv": scores["size_distribution"]["area_cv"],
            "boundary_regularity": scores["spatial_coherence"]["boundary_regularity"],
            "runtime_s": round(elapsed, 2),
        }

        if "detection" in scores:
            row["detection_f1"] = scores["detection"]["f1"]
            row["precision"] = scores["detection"]["precision"]
            row["recall"] = scores["detection"]["recall"]

        if "accuracy" in scores:
            row["mean_iou"] = scores["accuracy"]["mean_iou"]
            row["mean_dice"] = scores["accuracy"]["mean_dice"]
            row["ap_50"] = scores["accuracy"]["ap_50"]
            row["ap_50_95"] = scores["accuracy"]["ap_50_95"]

        if "marker_quality" in scores and len(scores["marker_quality"]) > 0:
            row["mean_snr"] = float(scores["marker_quality"]["snr"].mean())

        results.append(row)

        # Save mask if output_dir specified
        if output_dir is not None:
            _save_mask(mask, output_dir, method_name)

    return results


def run_benchmark(
    catalog: pd.DataFrame,
    config,
    output_dir: str | Path = "benchmark_results",
) -> pd.DataFrame:
    """Run the full benchmark across all ROIs in the catalog.

    Parameters
    ----------
    catalog
        DataFrame from ``build_benchmark_catalog()``.
    config
        ``BenchmarkConfig`` instance.
    output_dir
        Where to save results.

    Returns
    -------
    DataFrame with all results.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results_file = output_dir / "benchmark_results.csv"
    progress_file = output_dir / "progress.json"

    # Load previous progress for resume
    completed = set()
    if config.resume and progress_file.is_file():
        with open(progress_file) as f:
            completed = set(json.load(f).get("completed", []))
        logger.info("Resuming: %d ROIs already completed", len(completed))

    all_results = []

    # Load existing results
    if config.resume and results_file.is_file():
        existing = pd.read_csv(results_file)
        all_results = existing.to_dict("records")

    # Limit ROIs per dataset if configured
    if config.n_rois > 0:
        catalog = catalog.groupby("dataset").head(config.n_rois).reset_index(drop=True)

    total = len(catalog)
    for idx, roi_row in catalog.iterrows():
        roi_key = f"{roi_row['dataset']}_{roi_row['roi_id']}"

        if roi_key in completed:
            continue

        logger.info(
            "[%d/%d] Processing %s / %s",
            idx + 1,
            total,
            roi_row["dataset"],
            roi_row["roi_id"],
        )

        for strategy in config.strategies:
            methods = config.methods.get(strategy, [])
            if not methods:
                continue

            try:
                roi_results = run_single_roi(
                    tiff_path=roi_row["tiff_path"],
                    strategy=strategy,
                    methods=methods,
                    channel_csv_path=roi_row.get("channel_csv_path"),
                    prob_path=roi_row.get("prob_path"),
                    gt_mask_path=roi_row.get("mask_path"),
                    gpu=config.gpu,
                    prob_map_method=config.prob_map_method,
                    output_dir=output_dir / "masks" / roi_row["dataset"] / roi_row["roi_id"],
                )

                for r in roi_results:
                    r["dataset"] = roi_row["dataset"]
                    r["tissue"] = roi_row["tissue"]
                    r["sample_id"] = roi_row["sample_id"]
                    r["roi_id"] = roi_row["roi_id"]
                    r["has_gt"] = roi_row.get("has_gt", False)
                    all_results.append(r)

            except Exception as e:
                logger.error(
                    "Strategy %d failed for %s/%s: %s",
                    strategy,
                    roi_row["dataset"],
                    roi_row["roi_id"],
                    e,
                )

        # Mark as completed and save progress
        completed.add(roi_key)
        with open(progress_file, "w") as f:
            json.dump({"completed": list(completed)}, f)

        # Save intermediate results
        if len(all_results) > 0:
            pd.DataFrame(all_results).to_csv(results_file, index=False)

    df = pd.DataFrame(all_results) if all_results else pd.DataFrame()
    if len(df) > 0:
        df.to_csv(results_file, index=False)
        logger.info("Benchmark complete: %d results saved to %s", len(df), results_file)

    return df


def aggregate_results(results_df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Group results by strategy/method, dataset, and tissue type.

    Parameters
    ----------
    results_df
        Full benchmark results DataFrame.

    Returns
    -------
    Dict with keys:
    - ``"by_method"``: aggregated by (strategy, method)
    - ``"by_dataset"``: aggregated by (strategy, method, dataset)
    - ``"by_tissue"``: aggregated by (strategy, method, tissue)
    """
    if len(results_df) == 0:
        return {
            "by_method": pd.DataFrame(),
            "by_dataset": pd.DataFrame(),
            "by_tissue": pd.DataFrame(),
        }

    numeric_cols = results_df.select_dtypes(include=[np.number]).columns.tolist()
    agg_cols = [c for c in numeric_cols if c not in ("strategy",)]

    by_method = results_df.groupby(["strategy", "method"])[agg_cols].agg(["mean", "std", "count"])
    by_method.columns = ["_".join(col).strip() for col in by_method.columns]
    by_method = by_method.reset_index()

    by_dataset = results_df.groupby(["strategy", "method", "dataset"])[agg_cols].agg(
        ["mean", "std", "count"]
    )
    by_dataset.columns = ["_".join(col).strip() for col in by_dataset.columns]
    by_dataset = by_dataset.reset_index()

    by_tissue = results_df.groupby(["strategy", "method", "tissue"])[agg_cols].agg(
        ["mean", "std", "count"]
    )
    by_tissue.columns = ["_".join(col).strip() for col in by_tissue.columns]
    by_tissue = by_tissue.reset_index()

    return {
        "by_method": by_method,
        "by_dataset": by_dataset,
        "by_tissue": by_tissue,
    }


def _run_strategy(
    strategy: int,
    tiff_path: str | Path,
    methods: list[str],
    channel_csv_path: str | Path | None = None,
    panel_mapper=None,
    prob_path: str | Path | None = None,
    gpu: bool = False,
    prob_map_method: str = "gaussian",
) -> dict[str, np.ndarray]:
    """Dispatch to the appropriate strategy runner."""
    if strategy == 1:
        from sc_tools.bm.segment import run_all_strategy1

        # Load prob map or raw TIFF
        image = _load_image(prob_path if prob_path else tiff_path)
        return run_all_strategy1(image, methods=methods, gpu=gpu)

    elif strategy == 2:
        from sc_tools.bm.strategy_dna import run_all_strategy2

        return run_all_strategy2(
            tiff_path,
            methods=methods,
            channel_csv_path=channel_csv_path,
            panel_mapper=panel_mapper,
            prob_map_method=prob_map_method,
            gpu=gpu,
        )

    elif strategy == 3:
        from sc_tools.bm.strategy_hf import run_all_strategy3

        return run_all_strategy3(
            tiff_path,
            methods=methods,
            channel_csv_path=channel_csv_path,
            panel_mapper=panel_mapper,
            gpu=gpu,
        )

    elif strategy == 4:
        logger.warning("Strategy 4 (SegFormer) requires training first. Skipping in benchmark.")
        return {}

    else:
        raise ValueError(f"Unknown strategy: {strategy}")


def _load_image(path: str | Path) -> np.ndarray:
    """Load an image from a TIFF file."""
    import tifffile

    return tifffile.imread(str(path))


def _save_mask(mask: np.ndarray, output_dir: str | Path, method_name: str) -> None:
    """Save a mask to disk."""
    import tifffile

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"{method_name}_mask.tiff"
    tifffile.imwrite(str(path), mask.astype(np.int32))
