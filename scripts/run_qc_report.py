#!/usr/bin/env python3
"""
Single QC report script: produce all QC figures (2x2, 2x4, multipage, violin, scatter, HVG, SVG),
per-sample metrics with pass/fail classification, and a unified HTML report.

Usage (from repo root; sc_tools must be importable, e.g. pip install -e . or PYTHONPATH=.):
  python scripts/run_qc_report.py --raw projects/visium/ggo_visium/results/adata.annotated.p2.h5ad \\
    --post projects/visium/ggo_visium/results/adata.normalized.scored.p35.h5ad \\
    --figures-dir projects/visium/ggo_visium/figures \\
    --modality visium

From project dir (e.g. projects/visium/ggo_visium) with PYTHONPATH pointing at repo root:
  PYTHONPATH=/path/to/sc_tools python ../../../scripts/run_qc_report.py --raw results/adata.annotated.p2.h5ad \\
    --post results/adata.normalized.scored.p35.h5ad --figures-dir figures --modality visium

Snakemake (from repo root; uses Apptainer if run_apptainer.sh is used):
  snakemake -d projects/visium/ggo_visium -s projects/visium/ggo_visium/Snakefile figures/QC/qc_report.done

Output: figures/QC/raw/*.pdf|png, figures/QC/post/*.pdf|png, figures/QC/qc_report.html,
        figures/QC/sample_metrics.csv, figures/QC/qc_sample_pass.csv, figures/QC/qc_sample_fail.csv.
        Touches figures/QC/qc_report.done when done.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import scanpy as sc

from sc_tools.qc import (
    calculate_qc_metrics,
    classify_samples,
    compute_sample_metrics,
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
    figures_raw = figures_dir / "QC" / "raw"
    figures_post = figures_dir / "QC" / "post"
    figures_qc = figures_dir / "QC"
    figures_raw.mkdir(parents=True, exist_ok=True)
    figures_post.mkdir(parents=True, exist_ok=True)

    adata_raw = sc.read_h5ad(raw_path)
    adata_raw = adata_raw.copy()
    _ensure_qc_metrics(adata_raw)

    # ---- Sample-level metrics and classification ----
    sample_col = library_id_col
    if sample_col not in adata_raw.obs.columns:
        # Try common alternatives
        for alt in ["sample", "library_id", "sample_id"]:
            if alt in adata_raw.obs.columns:
                sample_col = alt
                break

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
        import matplotlib.pyplot as plt
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
    import matplotlib.pyplot as plt
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

    # ---- HTML report ----
    if metrics is not None and classified is not None:
        try:
            generate_qc_report(
                adata_raw, metrics, classified, figures_dir,
                output_path=figures_qc / "qc_report.html",
                sample_col=sample_col, modality=modality,
            )
        except ImportError:
            print("Warning: jinja2 not installed; skipping HTML report generation.", file=sys.stderr)

    if touch_done:
        done_file = figures_qc / "qc_report.done"
        done_file.parent.mkdir(parents=True, exist_ok=True)
        done_file.touch()


def main():
    parser = argparse.ArgumentParser(
        description="Produce QC report: figures, per-sample metrics, pass/fail, HTML report."
    )
    parser.add_argument("--raw", type=Path, required=True, help="Path to raw/pre-filter adata (e.g. p1 or p2).")
    parser.add_argument("--post", type=Path, default=None, help="Path to post-filter/normalized adata (e.g. p35). Optional.")
    parser.add_argument("--figures-dir", type=Path, default=Path("figures"), help="Base figures directory (default: figures).")
    parser.add_argument("--library-id-col", default="library_id", help="Obs column for library/sample (default: library_id).")
    parser.add_argument("--modality", default="visium", choices=["visium", "visium_hd", "xenium", "cosmx", "imc"],
                        help="Data modality (default: visium).")
    parser.add_argument("--dpi", type=int, default=300, help="DPI for PNG (default 300).")
    parser.add_argument("--mad-multiplier", type=float, default=None,
                        help="MAD multiplier for outlier detection (default: auto based on cohort size).")
    parser.add_argument("--thresholds", type=str, default=None,
                        help="JSON string or file path with threshold overrides.")
    parser.add_argument("--spaceranger-dirs", type=str, default=None,
                        help="JSON string or file path mapping sample -> Space Ranger outs/ dir.")
    parser.add_argument("--apply-filter", action="store_true",
                        help="Apply QC filter: write filtered p1.h5ad + backup.")
    parser.add_argument("--no-filter", action="store_true",
                        help="Skip filtering (report only). Overrides --apply-filter.")
    parser.add_argument("--results-dir", type=Path, default=None,
                        help="Results directory for filtered output (default: inferred from --raw).")
    parser.add_argument("--no-touch", action="store_true", help="Do not touch figures/QC/qc_report.done.")
    args = parser.parse_args()

    if not args.raw.exists():
        print(f"Error: raw adata not found: {args.raw}", file=sys.stderr)
        sys.exit(1)

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
        touch_done=not args.no_touch,
    )
    print(f"QC report written to {args.figures_dir / 'QC'}/")


if __name__ == "__main__":
    main()
