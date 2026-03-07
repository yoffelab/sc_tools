#!/usr/bin/env python
"""Run QC reports on ggo_visium, robin, and ggo-imc (panel_g + panel_h).

Prepares legacy ggo-imc files (adds library_id from roi, computes QC metrics)
and generates all applicable date-versioned QC reports for each project.
"""

import os
import sys
import traceback

import scanpy as sc

# Add repo root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sc_tools.qc import calculate_qc_metrics
from sc_tools.qc.report import (
    generate_post_filter_report,
    generate_post_integration_report,
    generate_pre_filter_report,
)
from sc_tools.qc.sample_qc import classify_samples, compute_sample_metrics

DATE_STAMP = "20260304"


def prepare_imc_adata(path, panel_name):
    """Load legacy IMC h5ad and add standard columns."""
    print(f"  Loading {panel_name}: {path}")
    adata = sc.read_h5ad(path)
    adata.obs_names_make_unique()

    # Add library_id from roi
    if "library_id" not in adata.obs.columns and "roi" in adata.obs.columns:
        adata.obs["library_id"] = adata.obs["roi"].astype(str)
        print(f"    Added library_id from roi ({adata.obs['library_id'].nunique()} unique)")

    # Add raw_data_dir from sample
    if "raw_data_dir" not in adata.obs.columns and "sample" in adata.obs.columns:
        adata.obs["raw_data_dir"] = adata.obs["sample"].astype(str)

    # Compute QC metrics if missing
    if "total_counts" not in adata.obs.columns:
        print("    Computing QC metrics...")
        calculate_qc_metrics(adata)
        print(f"    Done. total_counts range: {adata.obs['total_counts'].min():.0f}-{adata.obs['total_counts'].max():.0f}")

    return adata


def run_imc_reports(panel, raw_path, harmony_path, output_dir):
    """Run QC reports for one IMC panel."""
    print(f"\n{'='*60}")
    print(f"IMC ggo-imc {panel}")
    print(f"{'='*60}")

    # Pre-filter report (raw data)
    try:
        adata_raw = prepare_imc_adata(raw_path, f"{panel} raw")
        metrics = compute_sample_metrics(adata_raw, sample_col="roi")
        classified = classify_samples(metrics, modality="imc")

        print(f"  Generating pre-filter report...")
        out = generate_pre_filter_report(
            adata_raw, metrics, classified,
            output_dir=output_dir,
            modality="imc",
            sample_col="roi",
            date_stamp=DATE_STAMP,
        )
        print(f"  -> {out}")
    except Exception:
        print(f"  Pre-filter report FAILED:")
        traceback.print_exc()

    # Post-integration report (harmony data)
    try:
        adata_harmony = prepare_imc_adata(harmony_path, f"{panel} harmony")

        print(f"  Generating post-integration report...")
        out = generate_post_integration_report(
            adata_harmony,
            output_dir=output_dir,
            batch_key="roi",
            sample_col="roi",
            date_stamp=DATE_STAMP,
        )
        print(f"  -> {out}")
    except Exception:
        print(f"  Post-integration report FAILED:")
        traceback.print_exc()


def run_ggo_visium_reports():
    """Run QC reports for ggo_visium (no p1, use p2 as pre-filter equivalent)."""
    print(f"\n{'='*60}")
    print(f"Visium ggo_visium")
    print(f"{'='*60}")
    base = "projects/visium/ggo_visium"
    output_dir = f"{base}/figures/QC"
    os.makedirs(output_dir, exist_ok=True)

    # Use annotated as pre-filter (it's the earliest available checkpoint)
    _ann_new = f"{base}/results/adata.annotated.h5ad"
    _ann_old = f"{base}/results/adata.annotated.p2.h5ad"
    p2_path = _ann_new if os.path.exists(_ann_new) else _ann_old
    _norm_new = f"{base}/results/adata.normalized.h5ad"
    _norm_old = f"{base}/results/adata.normalized.p3.h5ad"
    p3_path = _norm_new if os.path.exists(_norm_new) else _norm_old

    # Pre-filter report using annotated (closest available)
    try:
        print(f"  Loading annotated: {p2_path}")
        adata_p2 = sc.read_h5ad(p2_path)
        adata_p2.obs_names_make_unique()

        # ggo_visium uses batch/library_id but not sample
        sample_col = "library_id" if "library_id" in adata_p2.obs.columns else "batch"
        if "raw_data_dir" not in adata_p2.obs.columns:
            adata_p2.obs["raw_data_dir"] = adata_p2.obs[sample_col].astype(str)

        metrics = compute_sample_metrics(adata_p2, sample_col=sample_col)
        classified = classify_samples(metrics, modality="visium")

        print(f"  Generating pre-filter report (from p2)...")
        out = generate_pre_filter_report(
            adata_p2, metrics, classified,
            output_dir=output_dir,
            modality="visium",
            sample_col=sample_col,
            date_stamp=DATE_STAMP,
        )
        print(f"  -> {out}")
    except Exception:
        print(f"  Pre-filter report FAILED:")
        traceback.print_exc()
        adata_p2 = None

    # Post-integration report
    try:
        print(f"  Loading p3: {p3_path}")
        adata_p3 = sc.read_h5ad(p3_path)
        adata_p3.obs_names_make_unique()

        batch_key = "library_id" if "library_id" in adata_p3.obs.columns else "batch"
        sample_col = batch_key

        print(f"  Generating post-integration report...")
        out = generate_post_integration_report(
            adata_p3,
            output_dir=output_dir,
            batch_key=batch_key,
            sample_col=sample_col,
            date_stamp=DATE_STAMP,
        )
        print(f"  -> {out}")
    except Exception:
        print(f"  Post-integration report FAILED:")
        traceback.print_exc()


def run_robin_reports():
    """Run all three QC reports for robin."""
    print(f"\n{'='*60}")
    print(f"Visium HD robin")
    print(f"{'='*60}")
    base = "projects/visium_hd/robin"
    output_dir = f"{base}/figures/QC"
    os.makedirs(output_dir, exist_ok=True)

    _raw_new = f"{base}/results/adata.raw.h5ad"
    _raw_old = f"{base}/results/adata.raw.p1.h5ad"
    p1_path = _raw_new if os.path.exists(_raw_new) else _raw_old
    _ann_new = f"{base}/results/adata.annotated.h5ad"
    _ann_old = f"{base}/results/adata.annotated.p2.h5ad"
    p2_path = _ann_new if os.path.exists(_ann_new) else _ann_old
    _norm_new = f"{base}/results/adata.normalized.h5ad"
    _norm_old = f"{base}/results/adata.normalized.p3.h5ad"
    p3_path = _norm_new if os.path.exists(_norm_new) else _norm_old

    adata_p1 = None

    # Pre-filter report
    try:
        print(f"  Loading raw: {p1_path}")
        adata_p1 = sc.read_h5ad(p1_path)
        adata_p1.obs_names_make_unique()

        sample_col = "sample" if "sample" in adata_p1.obs.columns else "batch"
        if "raw_data_dir" not in adata_p1.obs.columns:
            adata_p1.obs["raw_data_dir"] = adata_p1.obs[sample_col].astype(str)

        metrics = compute_sample_metrics(adata_p1, sample_col=sample_col)
        classified = classify_samples(metrics, modality="visium_hd")

        print(f"  Generating pre-filter report...")
        out = generate_pre_filter_report(
            adata_p1, metrics, classified,
            output_dir=output_dir,
            modality="visium_hd",
            sample_col=sample_col,
            date_stamp=DATE_STAMP,
        )
        print(f"  -> {out}")
    except Exception:
        print(f"  Pre-filter report FAILED:")
        traceback.print_exc()

    # Post-filter report (p1 vs p2)
    try:
        print(f"  Loading p2: {p2_path}")
        adata_p2 = sc.read_h5ad(p2_path)
        adata_p2.obs_names_make_unique()

        if adata_p1 is None:
            adata_p1 = sc.read_h5ad(p1_path)
            adata_p1.obs_names_make_unique()

        sample_col = "sample" if "sample" in adata_p2.obs.columns else "batch"
        metrics = compute_sample_metrics(adata_p2, sample_col=sample_col)
        classified = classify_samples(metrics, modality="visium_hd")

        print(f"  Generating post-filter report...")
        out = generate_post_filter_report(
            adata_p1, adata_p2, metrics, classified,
            output_dir=output_dir,
            modality="visium_hd",
            sample_col=sample_col,
            date_stamp=DATE_STAMP,
        )
        print(f"  -> {out}")
        del adata_p1, adata_p2  # free memory
    except Exception:
        print(f"  Post-filter report FAILED:")
        traceback.print_exc()

    # Post-integration report
    try:
        print(f"  Loading p3: {p3_path}")
        adata_p3 = sc.read_h5ad(p3_path)
        adata_p3.obs_names_make_unique()

        batch_key = "sample" if "sample" in adata_p3.obs.columns else "batch"

        print(f"  Generating post-integration report...")
        out = generate_post_integration_report(
            adata_p3,
            output_dir=output_dir,
            batch_key=batch_key,
            sample_col=batch_key,
            date_stamp=DATE_STAMP,
        )
        print(f"  -> {out}")
    except Exception:
        print(f"  Post-integration report FAILED:")
        traceback.print_exc()


if __name__ == "__main__":
    os.chdir("/Users/junbumkim/Documents/sc_tools")

    # 1. ggo-imc panel_g
    run_imc_reports(
        "panel_g",
        "projects/imc/ggo-imc/results/panel_g.raw.labeled.h5ad",
        "projects/imc/ggo-imc/results/panel_g.harmony.h5ad",
        "projects/imc/ggo-imc/figures/QC",
    )

    # 2. ggo-imc panel_h
    run_imc_reports(
        "panel_h",
        "projects/imc/ggo-imc/results/panel_h.raw.labeled.h5ad",
        "projects/imc/ggo-imc/results/panel_h.harmony.h5ad",
        "projects/imc/ggo-imc/figures/QC",
    )

    # 3. ggo_visium
    run_ggo_visium_reports()

    # 4. robin
    run_robin_reports()

    print(f"\n{'='*60}")
    print("All QC reports complete!")
    print(f"{'='*60}")
