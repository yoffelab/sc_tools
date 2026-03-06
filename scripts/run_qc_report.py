#!/usr/bin/env python3
"""
QC report script: produce QC figures and date-versioned HTML reports.

Supports three report types aligned with pipeline phases:
  - pre_filter:  Phase 1 entry (raw concatenated data)
  - post_filter: Phase 1-2 exit (after filtering + metadata)
  - post_integration: Phase 3 exit (after preprocessing)
  - all: generate all applicable reports

Usage:
  # Pre-filter only
  python scripts/run_qc_report.py --report pre_filter \\
    --adata results/adata.raw.p1.h5ad --figures-dir figures

  # Post-filter (pre + post)
  python scripts/run_qc_report.py --report post_filter \\
    --adata results/adata.raw.p1.h5ad --adata-post results/adata.annotated.p2.h5ad \\
    --figures-dir figures

  # Post-integration
  python scripts/run_qc_report.py --report post_integration \\
    --adata-integrated results/adata.normalized.p3.h5ad --figures-dir figures

  # Legacy (backward compat)
  python scripts/run_qc_report.py --raw results/adata.annotated.p2.h5ad \\
    --post results/adata.normalized.scored.p35.h5ad --figures-dir figures
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import scanpy as sc  # noqa: E402

from sc_tools.qc import (  # noqa: E402
    calculate_qc_metrics,
    classify_samples,
    compute_sample_metrics,
    generate_post_filter_report,
    generate_post_integration_report,
    generate_pre_filter_report,
    generate_qc_report,
    highly_variable_genes,
    plot_highly_variable_genes,
    plot_spatially_variable_genes,
    qc_2x2_grid,
    qc_2x4_pre_post,
    qc_sample_comparison_bar,
    qc_sample_scatter_matrix,
    qc_sample_violin_grouped,
    qc_scatter_counts_genes,
    qc_spatial_multipage,
    qc_violin_metrics,
    save_pass_fail_lists,
    spatially_variable_genes_per_library,
)


def _ensure_qc_metrics(adata):
    if "total_counts" in adata.obs.columns:
        return
    calculate_qc_metrics(adata, inplace=True, percent_top=(50, 100, 200, 500))


def _resolve_sample_col(adata, library_id_col):
    """Find a valid sample column in adata.obs."""
    if library_id_col in adata.obs.columns:
        return library_id_col
    for alt in ["sample", "library_id", "sample_id"]:
        if alt in adata.obs.columns:
            return alt
    return library_id_col


def _compute_metrics_classified(adata, sample_col, modality, mad_multiplier, thresholds, spaceranger_dirs):
    """Compute sample metrics and classification."""
    if sample_col not in adata.obs.columns:
        return None, None, sample_col
    metrics = compute_sample_metrics(
        adata, sample_col=sample_col, modality=modality,
        spaceranger_dirs=spaceranger_dirs,
    )
    classified = classify_samples(
        metrics, modality=modality, thresholds=thresholds,
        mad_multiplier=mad_multiplier if mad_multiplier is not None else 3.0,
    )
    return metrics, classified, sample_col


def run_pre_filter(
    adata_path, figures_dir, sample_col, modality, mad_multiplier, thresholds,
    spaceranger_dirs, segmentation_masks_dir, touch_done,
):
    """Generate pre-filter QC report."""
    import matplotlib.pyplot as plt

    adata = sc.read_h5ad(adata_path).copy()
    _ensure_qc_metrics(adata)
    sample_col = _resolve_sample_col(adata, sample_col)

    metrics, classified, sample_col = _compute_metrics_classified(
        adata, sample_col, modality, mad_multiplier, thresholds, spaceranger_dirs
    )

    if metrics is None or classified is None:
        print("Warning: could not compute sample metrics (sample column not found)", file=sys.stderr)
        return

    qc_dir = figures_dir / "QC"
    qc_dir.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(qc_dir / "sample_metrics.csv")
    save_pass_fail_lists(classified, qc_dir, sample_col=sample_col)

    # Also save standalone plots to figures/QC/raw/
    raw_dir = qc_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    qc_2x2_grid(adata, output_dir=raw_dir, basename="qc_2x2_raw", dpi=300)
    qc_violin_metrics(adata, output_dir=raw_dir, basename="qc_violin_raw", dpi=300)
    qc_scatter_counts_genes(adata, output_dir=raw_dir, basename="qc_scatter_raw", dpi=300)
    plt.close("all")

    generate_pre_filter_report(
        adata, metrics, classified, qc_dir,
        sample_col=sample_col, modality=modality,
        segmentation_masks_dir=segmentation_masks_dir,
    )

    if touch_done:
        (qc_dir / "pre_filter_qc.done").touch()


def run_post_filter(
    adata_pre_path, adata_post_path, figures_dir, sample_col, modality,
    mad_multiplier, thresholds, spaceranger_dirs, segmentation_masks_dir, touch_done,
):
    """Generate post-filter QC report."""
    import matplotlib.pyplot as plt

    adata_pre = sc.read_h5ad(adata_pre_path).copy()
    _ensure_qc_metrics(adata_pre)
    adata_post = sc.read_h5ad(adata_post_path).copy()
    _ensure_qc_metrics(adata_post)
    sample_col = _resolve_sample_col(adata_pre, sample_col)

    metrics, classified, sample_col = _compute_metrics_classified(
        adata_pre, sample_col, modality, mad_multiplier, thresholds, spaceranger_dirs
    )

    if metrics is None or classified is None:
        print("Warning: could not compute sample metrics", file=sys.stderr)
        return

    qc_dir = figures_dir / "QC"
    qc_dir.mkdir(parents=True, exist_ok=True)

    # Post standalone plots
    post_dir = qc_dir / "post"
    post_dir.mkdir(parents=True, exist_ok=True)
    qc_2x2_grid(adata_post, output_dir=post_dir, basename="qc_2x2_post", dpi=300)
    qc_violin_metrics(adata_post, output_dir=post_dir, basename="qc_violin_post", dpi=300)
    qc_2x4_pre_post(adata_pre, adata_post, output_dir=post_dir, basename="qc_2x4_pre_post", dpi=300)
    plt.close("all")

    generate_post_filter_report(
        adata_pre, adata_post, metrics, classified, qc_dir,
        sample_col=sample_col, modality=modality,
        segmentation_masks_dir=segmentation_masks_dir,
    )

    if touch_done:
        (qc_dir / "post_filter_qc.done").touch()


def run_post_integration(
    adata_path, figures_dir, sample_col, modality, batch_key, celltype_key,
    embedding_keys, segmentation_masks_dir, touch_done, run_benchmark=False,
):
    """Generate post-integration QC report."""
    adata = sc.read_h5ad(adata_path).copy()
    sample_col = _resolve_sample_col(adata, sample_col)

    qc_dir = figures_dir / "QC"
    qc_dir.mkdir(parents=True, exist_ok=True)

    # Optionally run full integration benchmark first
    if run_benchmark:
        try:
            from sc_tools.bm.integration import run_integration_benchmark

            adata, comparison_df = run_integration_benchmark(
                adata,
                modality=modality,
                batch_key=batch_key or "library_id",
                celltype_key=celltype_key,
            )
            print(f"Integration benchmark: {len(comparison_df)} methods evaluated")
            # Auto-detect embeddings from benchmark results
            if embedding_keys is None:
                from sc_tools.qc.report_utils import auto_detect_embeddings

                embedding_keys = auto_detect_embeddings(adata)
        except Exception as e:
            print(f"Warning: integration benchmark failed: {e}", file=sys.stderr)

    generate_post_integration_report(
        adata, qc_dir,
        embedding_keys=embedding_keys,
        batch_key=batch_key,
        celltype_key=celltype_key,
        sample_col=sample_col,
        modality=modality,
        segmentation_masks_dir=segmentation_masks_dir,
    )

    if touch_done:
        (qc_dir / "post_integration_qc.done").touch()


def run_segmentation(
    adata_path, masks_dir, figures_dir, sample_col, modality, touch_done,
):
    """Generate standalone segmentation QC report."""
    from sc_tools.qc import generate_segmentation_qc_report

    adata = sc.read_h5ad(adata_path).copy()
    sample_col = _resolve_sample_col(adata, sample_col)

    qc_dir = figures_dir / "QC"
    qc_dir.mkdir(parents=True, exist_ok=True)

    result = generate_segmentation_qc_report(
        adata, masks_dir, qc_dir,
        sample_col=sample_col, modality=modality,
    )

    if result is None:
        print(f"Warning: no masks found in {masks_dir}", file=sys.stderr)
        return

    if touch_done:
        (qc_dir / "segmentation_qc.done").touch()


def run_qc_report(
    raw_path: Path,
    post_path: Path | None,
    figures_dir: Path,
    library_id_col: str = "library_id",
    modality: str = "visium",
    dpi: int = 300,
    mad_multiplier: float | None = None,
    thresholds: dict | None = None,
    spaceranger_dirs: dict | None = None,
    apply_filter: bool = False,
    results_dir: Path | None = None,
    touch_done: bool = True,
) -> None:
    """Legacy run_qc_report function (backward compat)."""
    import matplotlib.pyplot as plt

    figures_raw = figures_dir / "QC" / "raw"
    figures_post = figures_dir / "QC" / "post"
    figures_qc = figures_dir / "QC"
    figures_raw.mkdir(parents=True, exist_ok=True)
    figures_post.mkdir(parents=True, exist_ok=True)

    adata_raw = sc.read_h5ad(raw_path)
    adata_raw = adata_raw.copy()
    _ensure_qc_metrics(adata_raw)

    # ---- Sample-level metrics and classification ----
    sample_col = _resolve_sample_col(adata_raw, library_id_col)

    metrics = None
    classified = None
    if sample_col in adata_raw.obs.columns:
        metrics = compute_sample_metrics(
            adata_raw,
            sample_col=sample_col,
            modality=modality,
            spaceranger_dirs=spaceranger_dirs,
        )
        classified = classify_samples(
            metrics,
            modality=modality,
            thresholds=thresholds,
            mad_multiplier=mad_multiplier if mad_multiplier is not None else 3.0,
        )
        # Save CSV outputs
        metrics.to_csv(figures_qc / "sample_metrics.csv")
        save_pass_fail_lists(classified, figures_qc, sample_col=sample_col)

        # Cross-sample comparison plots
        qc_sample_comparison_bar(
            metrics, classified=classified,
            output_dir=figures_qc, basename="qc_sample_comparison", dpi=dpi,
        )
        plt.close("all")

        qc_sample_violin_grouped(
            adata_raw, sample_col=sample_col, classified=classified,
            output_dir=figures_qc, basename="qc_sample_violin", dpi=dpi,
        )
        plt.close("all")

        qc_sample_scatter_matrix(
            metrics, classified=classified,
            output_dir=figures_qc, basename="qc_sample_scatter_matrix", dpi=dpi,
        )
        plt.close("all")

    # ---- Raw: 2x2, violin, scatter ----
    qc_2x2_grid(adata_raw, output_dir=figures_raw, basename="qc_2x2_raw", dpi=dpi)
    qc_violin_metrics(adata_raw, output_dir=figures_raw, basename="qc_violin_raw", dpi=dpi)
    qc_scatter_counts_genes(adata_raw, output_dir=figures_raw, basename="qc_scatter_raw", dpi=dpi)
    plt.close("all")

    # Multipage spatial (raw) when library_id and uns['spatial'] present
    if library_id_col in adata_raw.obs.columns and "spatial" in adata_raw.obsm:
        spatial_uns = adata_raw.uns.get("spatial", {})
        libs = adata_raw.obs[library_id_col].dropna().unique()
        if spatial_uns and all(str(lib) in spatial_uns for lib in libs):
            out_pdf = figures_raw / "qc_spatial_multipage_raw.pdf"
            qc_spatial_multipage(
                adata_raw, library_id_col, str(out_pdf), common_scale=True, dpi=dpi,
            )
            plt.close("all")

    # ---- Apply QC filter if requested ----
    if apply_filter and classified is not None and results_dir is not None:
        from sc_tools.qc import apply_qc_filter
        apply_qc_filter(
            adata_raw,
            classified,
            sample_col=sample_col,
            modality=modality,
            output_path=results_dir / "adata.raw.p1.h5ad",
            backup_path=results_dir / "adata.raw.p1.backup.qcfail_included.h5ad",
        )

    # ---- Post plots ----
    if post_path is not None and post_path.exists():
        adata_post = sc.read_h5ad(post_path)
        adata_post = adata_post.copy()
        _ensure_qc_metrics(adata_post)

        qc_2x2_grid(adata_post, output_dir=figures_post, basename="qc_2x2_post", dpi=dpi)
        qc_violin_metrics(adata_post, output_dir=figures_post, basename="qc_violin_post", dpi=dpi)
        qc_scatter_counts_genes(adata_post, output_dir=figures_post, basename="qc_scatter_post", dpi=dpi)
        plt.close("all")

        # 2x4 pre vs post
        qc_2x4_pre_post(
            adata_raw, adata_post,
            output_dir=figures_post, basename="qc_2x4_pre_post", dpi=dpi,
        )
        plt.close("all")

        # Multipage spatial (post)
        if library_id_col in adata_post.obs.columns and "spatial" in adata_post.obsm:
            spatial_uns = adata_post.uns.get("spatial", {})
            libs = adata_post.obs[library_id_col].dropna().unique()
            if spatial_uns and all(str(lib) in spatial_uns for lib in libs):
                out_pdf = figures_post / "qc_spatial_multipage_post.pdf"
                qc_spatial_multipage(
                    adata_post, library_id_col, str(out_pdf), common_scale=True, dpi=dpi,
                )
                plt.close("all")

        # HVG plot (post)
        adata_hvg = adata_post
        if "highly_variable" not in adata_post.var.columns:
            adata_hvg = adata_post.copy()
            sc.pp.normalize_total(adata_hvg, target_sum=1e4)
            sc.pp.log1p(adata_hvg)
            highly_variable_genes(adata_hvg, flavor="seurat", n_top_genes=2000, inplace=True)
        plot_highly_variable_genes(adata_hvg, output_dir=figures_post, basename="hvg_post", dpi=dpi)
        plt.close("all")

        # SVG: only when library_id present; skip without failing
        if library_id_col in adata_post.obs.columns and "spatial" in adata_post.obsm:
            per_lib = spatially_variable_genes_per_library(
                adata_post, library_id_col=library_id_col, n_top_genes=100,
            )
            if per_lib is not None:
                plot_spatially_variable_genes(
                    adata_post, output_dir=figures_post, basename="svg_post", dpi=dpi,
                )
                plt.close("all")

    # ---- HTML report (legacy) ----
    adata_post_for_report = None
    if post_path is not None and post_path.exists():
        try:
            adata_post_for_report = sc.read_h5ad(post_path)
            adata_post_for_report = adata_post_for_report.copy()
            _ensure_qc_metrics(adata_post_for_report)
        except Exception:
            adata_post_for_report = None

    if metrics is not None and classified is not None:
        try:
            generate_qc_report(
                adata_raw, metrics, classified, figures_dir,
                output_path=figures_qc / "qc_report.html",
                sample_col=sample_col, modality=modality,
                adata_post=adata_post_for_report,
            )
        except ImportError:
            print("Warning: jinja2 not installed; skipping HTML report generation.", file=sys.stderr)

    if touch_done:
        done_file = figures_qc / "qc_report.done"
        done_file.parent.mkdir(parents=True, exist_ok=True)
        done_file.touch()


def main():
    parser = argparse.ArgumentParser(
        description="Produce QC reports: date-versioned HTML with figures and metrics."
    )

    # New --report mode
    parser.add_argument(
        "--report", type=str, default=None,
        choices=["pre_filter", "post_filter", "post_integration", "segmentation", "all"],
        help="Report type (default: legacy mode if --raw provided, else 'all').",
    )

    # Input paths
    parser.add_argument("--adata", type=Path, default=None,
                        help="Path to AnnData for the report (pre-filter or post-integration).")
    parser.add_argument("--adata-post", type=Path, default=None,
                        help="Path to post-filter AnnData (for post_filter report).")
    parser.add_argument("--adata-integrated", type=Path, default=None,
                        help="Path to post-integration AnnData (for post_integration report).")

    # Legacy --raw/--post
    parser.add_argument("--raw", type=Path, default=None,
                        help="(Legacy) Path to raw/pre-filter adata.")
    parser.add_argument("--post", type=Path, default=None,
                        help="(Legacy) Path to post-filter/normalized adata.")

    # Common options
    parser.add_argument("--figures-dir", type=Path, default=Path("figures"),
                        help="Base figures directory (default: figures).")
    parser.add_argument("--library-id-col", default="library_id",
                        help="Obs column for library/sample (default: library_id).")
    parser.add_argument("--modality", default="visium",
                        choices=["visium", "visium_hd", "visium_hd_cell", "xenium", "cosmx", "imc"],
                        help="Data modality (default: visium).")
    parser.add_argument("--dpi", type=int, default=300, help="DPI for PNG (default 300).")
    parser.add_argument("--mad-multiplier", type=float, default=None,
                        help="MAD multiplier for outlier detection.")
    parser.add_argument("--thresholds", type=str, default=None,
                        help="JSON string or file path with threshold overrides.")
    parser.add_argument("--spaceranger-dirs", type=str, default=None,
                        help="JSON mapping sample -> Space Ranger outs/ dir.")

    # Post-integration specific
    parser.add_argument("--batch-key", type=str, default=None,
                        help="Batch column for integration metrics (auto-detect if absent).")
    parser.add_argument("--celltype-key", type=str, default=None,
                        help="Cell type column (optional; skip bio metrics if absent).")
    parser.add_argument("--embedding-keys", type=str, default=None,
                        help='JSON dict: {"name": "obsm_key"} (auto-detect if absent).')

    # Segmentation
    parser.add_argument("--segmentation-masks-dir", type=Path, default=None,
                        help="Directory with mask TIFFs for segmentation scoring.")

    # Integration benchmark
    parser.add_argument("--run-benchmark", action="store_true",
                        help="Run full integration benchmark (for post_integration report).")

    # Filter options (legacy)
    parser.add_argument("--apply-filter", action="store_true",
                        help="Apply QC filter: write filtered p1.h5ad + backup.")
    parser.add_argument("--no-filter", action="store_true",
                        help="Skip filtering. Overrides --apply-filter.")
    parser.add_argument("--results-dir", type=Path, default=None,
                        help="Results directory for filtered output.")
    parser.add_argument("--no-touch", action="store_true",
                        help="Do not touch sentinel done files.")

    args = parser.parse_args()

    # Parse JSON args
    thresholds = None
    if args.thresholds:
        p = Path(args.thresholds)
        if p.exists():
            thresholds = json.loads(p.read_text())
        else:
            thresholds = json.loads(args.thresholds)

    spaceranger_dirs = None
    if args.spaceranger_dirs:
        p = Path(args.spaceranger_dirs)
        if p.exists():
            spaceranger_dirs = json.loads(p.read_text())
        else:
            spaceranger_dirs = json.loads(args.spaceranger_dirs)

    embedding_keys = None
    if args.embedding_keys:
        embedding_keys = json.loads(args.embedding_keys)

    touch_done = not args.no_touch
    seg_dir = str(args.segmentation_masks_dir) if args.segmentation_masks_dir else None

    # Legacy mode: --raw provided without --report
    if args.report is None and args.raw is not None:
        if not args.raw.exists():
            print(f"Error: raw adata not found: {args.raw}", file=sys.stderr)
            sys.exit(1)
        do_filter = args.apply_filter and not args.no_filter
        results_dir = args.results_dir
        if results_dir is None and do_filter:
            results_dir = args.raw.parent
        run_qc_report(
            raw_path=args.raw,
            post_path=args.post,
            figures_dir=args.figures_dir,
            library_id_col=args.library_id_col,
            modality=args.modality,
            dpi=args.dpi,
            mad_multiplier=args.mad_multiplier,
            thresholds=thresholds,
            spaceranger_dirs=spaceranger_dirs,
            apply_filter=do_filter,
            results_dir=results_dir,
            touch_done=touch_done,
        )
        print(f"QC report written to {args.figures_dir / 'QC'}/")
        return

    # New report mode
    report_type = args.report or "all"

    if report_type in ("pre_filter", "all"):
        adata_path = args.adata or args.raw
        if adata_path is None:
            print("Error: --adata (or --raw) required for pre_filter report", file=sys.stderr)
            sys.exit(1)
        if not adata_path.exists():
            print(f"Error: adata not found: {adata_path}", file=sys.stderr)
            sys.exit(1)
        run_pre_filter(
            adata_path, args.figures_dir, args.library_id_col, args.modality,
            args.mad_multiplier, thresholds, spaceranger_dirs, seg_dir, touch_done,
        )
        print(f"Pre-filter QC report written to {args.figures_dir / 'QC'}/")

    if report_type in ("post_filter", "all"):
        adata_pre_path = args.adata or args.raw
        adata_post_path = args.adata_post or args.post
        if adata_pre_path is None or adata_post_path is None:
            if report_type == "post_filter":
                print("Error: --adata and --adata-post required for post_filter report", file=sys.stderr)
                sys.exit(1)
            else:
                print("Skipping post_filter report (--adata-post not provided)", file=sys.stderr)
        elif adata_pre_path.exists() and adata_post_path.exists():
            run_post_filter(
                adata_pre_path, adata_post_path, args.figures_dir, args.library_id_col,
                args.modality, args.mad_multiplier, thresholds, spaceranger_dirs,
                seg_dir, touch_done,
            )
            print(f"Post-filter QC report written to {args.figures_dir / 'QC'}/")

    if report_type in ("post_integration", "all"):
        adata_int_path = args.adata_integrated or args.adata
        if adata_int_path is None:
            if report_type == "post_integration":
                print("Error: --adata-integrated (or --adata) required", file=sys.stderr)
                sys.exit(1)
            else:
                print("Skipping post_integration report (--adata-integrated not provided)", file=sys.stderr)
        elif adata_int_path.exists():
            run_post_integration(
                adata_int_path, args.figures_dir, args.library_id_col, args.modality,
                args.batch_key, args.celltype_key, embedding_keys, seg_dir, touch_done,
                run_benchmark=args.run_benchmark,
            )
            print(f"Post-integration QC report written to {args.figures_dir / 'QC'}/")

    if report_type in ("segmentation", "all"):
        adata_path = args.adata or args.raw
        if args.segmentation_masks_dir is None:
            if report_type == "segmentation":
                print("Error: --segmentation-masks-dir required for segmentation report", file=sys.stderr)
                sys.exit(1)
            # Skip silently in "all" mode
        elif adata_path is None:
            if report_type == "segmentation":
                print("Error: --adata required for segmentation report", file=sys.stderr)
                sys.exit(1)
        elif adata_path.exists():
            run_segmentation(
                adata_path, args.segmentation_masks_dir, args.figures_dir,
                args.library_id_col, args.modality, touch_done,
            )
            print(f"Segmentation QC report written to {args.figures_dir / 'QC'}/")

    # Legacy umbrella sentinel
    if report_type == "all" and touch_done:
        qc_dir = args.figures_dir / "QC"
        qc_dir.mkdir(parents=True, exist_ok=True)
        (qc_dir / "qc_report.done").touch()


if __name__ == "__main__":
    main()
