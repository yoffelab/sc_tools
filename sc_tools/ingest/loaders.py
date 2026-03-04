"""Modality-specific AnnData loaders and sample concatenation.

Provides standardized loading functions for Visium, Visium HD, Xenium,
and IMC data. Each loader sets required Phase 1 obs/obsm keys.
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def load_visium_sample(
    spaceranger_dir: str | Path,
    sample_id: str,
    *,
    load_images: bool = True,
) -> ad.AnnData:
    """Load one Visium sample from Space Ranger output.

    Sets obs['sample'], obs['library_id'], obs['raw_data_dir'],
    and obsm['spatial'].

    Parameters
    ----------
    spaceranger_dir
        Path to Space Ranger output directory (containing ``outs/``).
    sample_id
        Sample identifier to store in obs['sample'].
    load_images
        Whether to load H&E images into the AnnData object.

    Returns
    -------
    AnnData with spatial coordinates and sample metadata.
    """
    import scanpy as sc

    spaceranger_dir = Path(spaceranger_dir)
    outs_dir = spaceranger_dir / "outs" if (spaceranger_dir / "outs").exists() else spaceranger_dir

    adata = sc.read_visium(outs_dir, load_images=load_images)
    adata.obs["sample"] = sample_id
    adata.obs["library_id"] = sample_id
    adata.obs["raw_data_dir"] = str(spaceranger_dir)
    adata.var_names_make_unique()

    logger.info(
        "Loaded Visium sample %s: %d spots x %d genes",
        sample_id,
        adata.n_obs,
        adata.n_vars,
    )
    return adata


def load_visium_hd_sample(
    spaceranger_dir: str | Path,
    sample_id: str,
    *,
    bin_size: str = "square_008um",
    load_images: bool = False,
) -> ad.AnnData:
    """Load one Visium HD sample from binned Space Ranger output.

    Reads tissue positions from parquet and sets spatial coordinates.

    Parameters
    ----------
    spaceranger_dir
        Path to Space Ranger output directory containing binned outputs.
    sample_id
        Sample identifier.
    bin_size
        Bin size subdirectory name (default: square_008um).
    load_images
        Whether to load images.

    Returns
    -------
    AnnData with spatial coordinates and sample metadata.
    """
    try:
        import squidpy as sq
    except ImportError as e:
        raise ImportError("squidpy is required for Visium HD loading: pip install squidpy") from e

    spaceranger_dir = Path(spaceranger_dir)
    bin_dir = spaceranger_dir / bin_size

    if not bin_dir.exists():
        # Try under outs/
        bin_dir = spaceranger_dir / "outs" / bin_size
    if not bin_dir.exists():
        raise FileNotFoundError(
            f"Bin directory not found: {spaceranger_dir / bin_size} "
            f"or {spaceranger_dir / 'outs' / bin_size}"
        )

    adata = sq.read.visium(
        bin_dir,
        load_images=load_images,
        counts_file="filtered_feature_bc_matrix.h5",
    )

    # Load tissue positions from parquet
    pos_file = bin_dir / "spatial" / "tissue_positions.parquet"
    if pos_file.exists():
        df_pos = pd.read_parquet(pos_file).set_index("barcode")
        # Keep only barcodes present in adata
        common = adata.obs_names.intersection(df_pos.index)
        if len(common) > 0:
            adata = adata[common].copy()
            adata.obs = df_pos.loc[common]

    # Set spatial coordinates
    if "pxl_col_in_fullres" in adata.obs.columns:
        adata.obsm["spatial"] = np.array(adata.obs[["pxl_col_in_fullres", "pxl_row_in_fullres"]])

    adata.obs["sample"] = sample_id
    adata.obs["library_id"] = sample_id
    adata.obs["raw_data_dir"] = str(spaceranger_dir)
    adata.var_names_make_unique()

    logger.info(
        "Loaded Visium HD sample %s (%s): %d spots x %d genes",
        sample_id,
        bin_size,
        adata.n_obs,
        adata.n_vars,
    )
    return adata


def load_xenium_sample(
    xenium_dir: str | Path,
    sample_id: str,
) -> ad.AnnData:
    """Load one Xenium sample from output directory.

    Attempts spatialdata-io first, falls back to scanpy read_10x_h5.

    Parameters
    ----------
    xenium_dir
        Path to Xenium output directory.
    sample_id
        Sample identifier.

    Returns
    -------
    AnnData with spatial coordinates.
    """
    xenium_dir = Path(xenium_dir)

    try:
        import spatialdata_io

        sdata = spatialdata_io.xenium(xenium_dir)
        # Extract the table (AnnData) from SpatialData
        if hasattr(sdata, "tables") and "table" in sdata.tables:
            adata = sdata.tables["table"]
        elif hasattr(sdata, "table"):
            adata = sdata.table
        else:
            raise ValueError("Could not extract AnnData from SpatialData object")
    except ImportError:
        import scanpy as sc

        h5_path = xenium_dir / "cell_feature_matrix.h5"
        if not h5_path.exists():
            h5_path = xenium_dir / "cell_feature_matrix" / "matrix.mtx.gz"
        adata = (
            sc.read_10x_h5(str(h5_path))
            if h5_path.exists()
            else sc.read_10x_mtx(str(xenium_dir / "cell_feature_matrix"))
        )

        # Load coordinates
        coords_file = xenium_dir / "cells.csv.gz"
        if not coords_file.exists():
            coords_file = xenium_dir / "cells.csv"
        if coords_file.exists():
            coords = pd.read_csv(coords_file)
            if "x_centroid" in coords.columns:
                adata.obsm["spatial"] = coords[["x_centroid", "y_centroid"]].values

    adata.obs["sample"] = sample_id
    adata.obs["library_id"] = sample_id
    adata.obs["raw_data_dir"] = str(xenium_dir)
    adata.var_names_make_unique()

    logger.info(
        "Loaded Xenium sample %s: %d cells x %d genes",
        sample_id,
        adata.n_obs,
        adata.n_vars,
    )
    return adata


def load_imc_sample(
    h5ad_path: str | Path,
    sample_id: str,
) -> ad.AnnData:
    """Load one IMC ROI from an h5ad file.

    Parameters
    ----------
    h5ad_path
        Path to .h5ad file for a single ROI.
    sample_id
        Sample/ROI identifier.

    Returns
    -------
    AnnData with sample annotation.
    """
    adata = ad.read_h5ad(h5ad_path)
    adata.obs["sample"] = sample_id
    adata.obs["library_id"] = sample_id
    adata.obs["raw_data_dir"] = str(h5ad_path)
    adata.var_names_make_unique()

    # Ensure spatial coordinates exist
    if "spatial" not in adata.obsm and "X_spatial" in adata.obsm:
        adata.obsm["spatial"] = adata.obsm["X_spatial"]

    logger.info(
        "Loaded IMC sample %s: %d cells x %d markers",
        sample_id,
        adata.n_obs,
        adata.n_vars,
    )
    return adata


def concat_samples(
    adatas: list[ad.AnnData],
    *,
    sample_col: str = "sample",
    calculate_qc: bool = True,
) -> ad.AnnData:
    """Concatenate multiple sample AnnDatas with proper handling.

    Parameters
    ----------
    adatas
        List of AnnData objects to concatenate.
    sample_col
        Column in obs used as sample identifier.
    calculate_qc
        If True, run scanpy calculate_qc_metrics after concatenation.

    Returns
    -------
    Concatenated AnnData with QC metrics.
    """
    if not adatas:
        raise ValueError("No AnnData objects to concatenate")

    if len(adatas) == 1:
        adata = adatas[0].copy()
    else:
        adata = ad.concat(adatas, merge="same", uns_merge="same")

    # Ensure spatial coordinates are preserved as numpy array
    if "spatial" in adata.obsm:
        adata.obsm["spatial"] = np.array(adata.obsm["spatial"])

    if calculate_qc:
        import scanpy as sc

        sc.pp.calculate_qc_metrics(adata, inplace=True)

    logger.info(
        "Concatenated %d samples: %d total obs x %d vars",
        len(adatas),
        adata.n_obs,
        adata.n_vars,
    )
    return adata
