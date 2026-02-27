#!/usr/bin/env python3
"""
Single QC report script: produce all QC figures (2x2, 2x4, multipage, violin, scatter, HVG, SVG)
without writing any AnnData. Run from repo root or project dir.

Usage (from repo root; sc_tools must be importable, e.g. pip install -e . or PYTHONPATH=.):
  python scripts/run_qc_report.py --raw projects/visium/ggo_visium/results/adata.annotated.p2.h5ad \\
    --post projects/visium/ggo_visium/results/adata.normalized.scored.p35.h5ad \\
    --figures-dir projects/visium/ggo_visium/figures

From project dir (e.g. projects/visium/ggo_visium) with PYTHONPATH pointing at repo root:
  PYTHONPATH=/path/to/sc_tools python ../../../scripts/run_qc_report.py --raw results/adata.annotated.p2.h5ad \\
    --post results/adata.normalized.scored.p35.h5ad --figures-dir figures

Snakemake (from repo root; uses Apptainer if run_apptainer.sh is used):
  snakemake -d projects/visium/ggo_visium -s projects/visium/ggo_visium/Snakefile figures/QC/qc_report.done

Output: figures/QC/raw/*.pdf|png and figures/QC/post/*.pdf|png. Touches figures/QC/qc_report.done when done.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import scanpy as sc

from sc_tools.qc import (
    calculate_qc_metrics,
    highly_variable_genes,
    plot_highly_variable_genes,
    plot_spatially_variable_genes,
    qc_2x2_grid,
    qc_2x4_pre_post,
    qc_scatter_counts_genes,
    qc_spatial_multipage,
    qc_violin_metrics,
    spatially_variable_genes_per_library,
)


def _ensure_qc_metrics(adata):
    if "total_counts" in adata.obs.columns:
        return
    calculate_qc_metrics(adata, inplace=True, percent_top=(50, 100, 200, 500))


def run_qc_report(
    raw_path: Path,
    post_path: Path,
    figures_dir: Path,
    library_id_col: str = "library_id",
    dpi: int = 300,
    touch_done: bool = True,
) -> None:
    figures_raw = figures_dir / "QC" / "raw"
    figures_post = figures_dir / "QC" / "post"
    figures_raw.mkdir(parents=True, exist_ok=True)
    figures_post.mkdir(parents=True, exist_ok=True)

    adata_raw = sc.read_h5ad(raw_path)
    adata_raw = adata_raw.copy()
    _ensure_qc_metrics(adata_raw)

    adata_post = sc.read_h5ad(post_path)
    adata_post = adata_post.copy()
    _ensure_qc_metrics(adata_post)

    # Raw: 2x2, violin, scatter
    qc_2x2_grid(adata_raw, output_dir=figures_raw, basename="qc_2x2_raw", dpi=dpi)
    qc_violin_metrics(adata_raw, output_dir=figures_raw, basename="qc_violin_raw", dpi=dpi)
    qc_scatter_counts_genes(adata_raw, output_dir=figures_raw, basename="qc_scatter_raw", dpi=dpi)

    # Post: 2x2, violin, scatter
    qc_2x2_grid(adata_post, output_dir=figures_post, basename="qc_2x2_post", dpi=dpi)
    qc_violin_metrics(adata_post, output_dir=figures_post, basename="qc_violin_post", dpi=dpi)
    qc_scatter_counts_genes(adata_post, output_dir=figures_post, basename="qc_scatter_post", dpi=dpi)

    # 2x4 pre vs post
    qc_2x4_pre_post(
        adata_raw,
        adata_post,
        output_dir=figures_post,
        basename="qc_2x4_pre_post",
        dpi=dpi,
    )

    # Multipage spatial (raw) when library_id and uns['spatial'] present
    if library_id_col in adata_raw.obs.columns and "spatial" in adata_raw.obsm:
        spatial_uns = adata_raw.uns.get("spatial", {})
        libs = adata_raw.obs[library_id_col].dropna().unique()
        if spatial_uns and all(str(l) in spatial_uns for l in libs):
            out_pdf = figures_raw / "qc_spatial_multipage_raw.pdf"
            qc_spatial_multipage(
                adata_raw,
                library_id_col,
                str(out_pdf),
                common_scale=True,
                dpi=dpi,
            )

    # Multipage spatial (post)
    if library_id_col in adata_post.obs.columns and "spatial" in adata_post.obsm:
        spatial_uns = adata_post.uns.get("spatial", {})
        libs = adata_post.obs[library_id_col].dropna().unique()
        if spatial_uns and all(str(l) in spatial_uns for l in libs):
            out_pdf = figures_post / "qc_spatial_multipage_post.pdf"
            qc_spatial_multipage(
                adata_post,
                library_id_col,
                str(out_pdf),
                common_scale=True,
                dpi=dpi,
            )

    # HVG plot (post): run HVG on copy if not present (seurat flavor avoids skmisc dependency)
    adata_hvg = adata_post
    if "highly_variable" not in adata_post.var.columns:
        adata_hvg = adata_post.copy()
        sc.pp.normalize_total(adata_hvg, target_sum=1e4)
        sc.pp.log1p(adata_hvg)
        highly_variable_genes(adata_hvg, flavor="seurat", n_top_genes=2000, inplace=True)
    plot_highly_variable_genes(
        adata_hvg,
        output_dir=figures_post,
        basename="hvg_post",
        dpi=dpi,
    )

    # SVG: only when library_id present; skip without failing
    if library_id_col in adata_post.obs.columns and "spatial" in adata_post.obsm:
        per_lib = spatially_variable_genes_per_library(
            adata_post,
            library_id_col=library_id_col,
            n_top_genes=100,
        )
        if per_lib is not None:
            plot_spatially_variable_genes(
                adata_post,
                output_dir=figures_post,
                basename="svg_post",
                dpi=dpi,
            )

    if touch_done:
        done_file = figures_dir / "QC" / "qc_report.done"
        done_file.parent.mkdir(parents=True, exist_ok=True)
        done_file.touch()


def main():
    parser = argparse.ArgumentParser(description="Produce QC report figures (no .h5ad written).")
    parser.add_argument("--raw", type=Path, required=True, help="Path to raw/pre-filter adata (e.g. p1 or p2).")
    parser.add_argument("--post", type=Path, required=True, help="Path to post-filter/normalized adata (e.g. p35).")
    parser.add_argument("--figures-dir", type=Path, default=Path("figures"), help="Base figures directory (default: figures).")
    parser.add_argument("--library-id-col", default="library_id", help="Obs column for library/sample (default: library_id).")
    parser.add_argument("--dpi", type=int, default=300, help="DPI for PNG (default 300).")
    parser.add_argument("--no-touch", action="store_true", help="Do not touch figures/QC/qc_report.done.")
    args = parser.parse_args()

    if not args.raw.exists():
        print(f"Error: raw adata not found: {args.raw}", file=sys.stderr)
        sys.exit(1)
    if not args.post.exists():
        print(f"Error: post adata not found: {args.post}", file=sys.stderr)
        sys.exit(1)

    run_qc_report(
        raw_path=args.raw,
        post_path=args.post,
        figures_dir=args.figures_dir,
        library_id_col=args.library_id_col,
        dpi=args.dpi,
        touch_done=not args.no_touch,
    )
    print(f"QC report written to {args.figures_dir / 'QC'}/raw/ and .../post/")


if __name__ == "__main__":
    main()
