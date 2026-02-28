"""
Multipage spatial plot of cell-type proportions from deconvolution.

One page per library_id, one subplot per cell type showing the proportion
values overlaid on spatial coordinates.

Usage (from project root):
    python <repo_root>/scripts/plot_deconvolution_spatial.py \
        --input results/adata.deconvolution.h5ad \
        --output figures/deconvolution/spatial_proportions.pdf \
        --batch-key library_id
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.backends.backend_pdf import PdfPages

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


def plot_deconvolution_spatial(
    adata_path: str,
    output_path: str,
    batch_key: str = "library_id",
    max_celltypes_per_page: int = 12,
    spot_size: float | None = None,
    figsize_per_panel: tuple[float, float] = (4, 4),
    dpi: int = 200,
) -> None:
    """Generate multipage spatial PDF of cell-type proportions.

    Parameters
    ----------
    adata_path
        Path to AnnData with ``obsm['cell_type_proportions']``.
    output_path
        PDF output path.
    batch_key
        Column in obs for library/sample grouping.
    max_celltypes_per_page
        Max cell types per page (grid wraps).
    spot_size
        Spot size for scatter. None = auto.
    figsize_per_panel
        Size per subplot panel.
    dpi
        Output DPI.
    """
    logger.info(f"Loading {adata_path}")
    adata = sc.read_h5ad(adata_path)

    if "cell_type_proportions" not in adata.obsm:
        raise ValueError("obsm['cell_type_proportions'] not found in AnnData")

    props = adata.obsm["cell_type_proportions"]
    if isinstance(props, pd.DataFrame):
        ct_names = list(props.columns)
        props_arr = props.values
    else:
        props_arr = np.asarray(props)
        ct_names = [f"celltype_{i}" for i in range(props_arr.shape[1])]

    logger.info(f"Cell types: {len(ct_names)}")
    logger.info(f"Cell type names: {ct_names}")

    if batch_key not in adata.obs.columns:
        logger.info(f"batch_key {batch_key!r} not found, treating as single batch")
        adata.obs["_batch"] = "all"
        batch_key = "_batch"

    library_ids = sorted(adata.obs[batch_key].unique())
    logger.info(f"Libraries: {len(library_ids)}: {library_ids}")

    # Check if spatial coordinates are available
    has_spatial = "spatial" in adata.obsm
    if not has_spatial:
        logger.warning("No spatial coordinates found in obsm['spatial'], using UMAP or random layout")

    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # PNG output directory (same folder as PDF)
    png_dir = out_path.parent
    png_dpi = 300

    with PdfPages(str(out_path)) as pdf:
        for lib_id in library_ids:
            lib_mask = adata.obs[batch_key] == lib_id
            lib_idx = np.where(lib_mask)[0]
            n_spots = lib_mask.sum()

            if n_spots == 0:
                continue

            # Get coordinates
            if has_spatial:
                coords = adata.obsm["spatial"][lib_idx]
                x, y = coords[:, 0], coords[:, 1]
            elif "X_umap" in adata.obsm:
                coords = adata.obsm["X_umap"][lib_idx]
                x, y = coords[:, 0], coords[:, 1]
            else:
                x = np.random.default_rng(42).uniform(0, 1, n_spots)
                y = np.random.default_rng(42).uniform(0, 1, n_spots)

            lib_props = props_arr[lib_idx]

            # Layout: grid of subplots
            n_cts = len(ct_names)
            ncols = min(4, n_cts)
            nrows = int(np.ceil(n_cts / ncols))
            fig_w = figsize_per_panel[0] * ncols
            fig_h = figsize_per_panel[1] * nrows

            fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h))
            if nrows == 1 and ncols == 1:
                axes = np.array([[axes]])
            elif nrows == 1:
                axes = axes[np.newaxis, :]
            elif ncols == 1:
                axes = axes[:, np.newaxis]

            fig.suptitle(f"Cell-type proportions: {lib_id} ({n_spots} spots)", fontsize=14, y=1.01)

            # Auto spot size: scale inversely with spot count
            # Visium (~3500 spots) -> ~34, Visium HD (~50K) -> ~2.4, HD (~150K) -> ~0.8
            s = spot_size if spot_size is not None else max(0.3, min(50, 120000 / n_spots))

            for idx, ct in enumerate(ct_names):
                row, col = divmod(idx, ncols)
                ax = axes[row, col]
                ax.set_facecolor("white")
                values = lib_props[:, idx]

                # Use 98th percentile for vmax to improve contrast
                # (outlier high-proportion spots otherwise wash out the rest)
                p98 = np.percentile(values, 98)
                vmax = max(p98, 0.01)

                scatter = ax.scatter(
                    x, y, c=values, s=s, cmap="viridis",
                    vmin=0, vmax=vmax,
                    alpha=0.8, edgecolors="none",
                )
                ax.set_title(ct, fontsize=9)
                ax.set_aspect("equal")
                ax.invert_yaxis()
                ax.set_xticks([])
                ax.set_yticks([])
                plt.colorbar(scatter, ax=ax, shrink=0.6, pad=0.02)

            # Hide unused axes
            for idx in range(n_cts, nrows * ncols):
                row, col = divmod(idx, ncols)
                axes[row, col].set_visible(False)

            plt.tight_layout()
            pdf.savefig(fig, dpi=dpi, bbox_inches="tight")

            # Save per-library PNG at 300 DPI
            png_path = png_dir / f"{lib_id}.png"
            fig.savefig(png_path, dpi=png_dpi, bbox_inches="tight")
            logger.info(f"  Page: {lib_id} ({n_spots} spots, {n_cts} cell types) -> {png_path}")

            plt.close(fig)

    logger.info(f"Saved PDF: {out_path}")
    logger.info(f"Saved PNGs: {png_dir}/*.png ({len(library_ids)} files, {png_dpi} DPI)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Spatial plot of deconvolution proportions")
    parser.add_argument("--input", required=True, help="Path to adata.deconvolution.h5ad")
    parser.add_argument("--output", required=True, help="Output PDF path")
    parser.add_argument("--batch-key", default="library_id", help="Batch/library column in obs")
    parser.add_argument("--spot-size", type=float, default=None, help="Spot size for scatter")
    parser.add_argument("--dpi", type=int, default=200, help="Output DPI")
    args = parser.parse_args()

    plot_deconvolution_spatial(
        adata_path=args.input,
        output_path=args.output,
        batch_key=args.batch_key,
        spot_size=args.spot_size,
        dpi=args.dpi,
    )
