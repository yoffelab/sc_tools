"""Modality-specific AnnData loaders and sample concatenation (Phase 0b).

Provides standardized loading functions for Visium, Visium HD, Xenium,
IMC, and CosMx data. Each loader reads from the Phase 0a platform output
directory and writes a per-sample AnnData with standardized obs/obsm keys
(sample, library_id, raw_data_dir, spatial).
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Internal helper: remote URI -> local path
# ---------------------------------------------------------------------------


def _local_path(uri: str | os.PathLike) -> Path:
    """Resolve *uri* to a local path (downloads remote files to tmp).

    For local paths this is a no-op.  For remote URIs (s3://, sftp://,
    etc.) the file is downloaded to a temporary location.  The caller is
    responsible for using the path before the temporary file is cleaned up
    (i.e. use only within a ``with_local_copy()`` context or read
    immediately).

    Note: directory-based loaders (Visium, Xenium) require local directories
    and do not benefit from URI resolution here.  For those, pass a local
    path.  Only IMC ``cells.h5ad`` files are fetched remotely.
    """
    try:
        from sc_tools.storage import _is_local, resolve_fs

        fs, path = resolve_fs(str(uri))
        if _is_local(fs):
            return Path(path)
        # Remote: not supported for directory-based loaders
        raise ValueError(
            f"Remote URIs are not supported for directory-based loaders. "
            f"Use sc_tools.storage.with_local_copy() to download first, "
            f"then pass the local path. URI: {uri}"
        )
    except ImportError:
        return Path(str(uri))


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

    # Search common bin directory locations
    bin_candidates = [
        spaceranger_dir / bin_size,
        spaceranger_dir / "outs" / bin_size,
        spaceranger_dir / "outs" / "binned_outputs" / bin_size,
    ]
    bin_dir = next((p for p in bin_candidates if p.exists()), None)
    if bin_dir is None:
        raise FileNotFoundError(
            f"Bin directory not found for {bin_size}. Tried: {[str(p) for p in bin_candidates]}"
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


def _extract_centroids_from_geojson(
    geojson_path: Path,
    obs_names: pd.Index | None = None,
) -> pd.DataFrame:
    """Extract cell centroids from a SpaceRanger 4 cell_segmentations.geojson.

    Each GeoJSON Feature has a Polygon geometry and a ``cell_id`` property
    (integer). The centroid is computed as the mean of the polygon exterior
    coordinates.

    The index format is auto-detected to match ``obs_names``. SpaceRanger 4
    uses ``cellid_XXXXXXXXX-1`` (zero-padded 9-digit cell_id with ``-1``
    suffix). If ``obs_names`` is provided and its first element matches this
    pattern, the index is formatted accordingly; otherwise plain string
    cell_id is used.

    Returns a DataFrame with columns ``x_centroid``, ``y_centroid``.
    """
    import json

    with open(geojson_path) as fh:
        data = json.load(fh)

    # Detect index format from obs_names
    use_cellid_format = False
    if obs_names is not None and len(obs_names) > 0:
        first = str(obs_names[0])
        if first.startswith("cellid_"):
            use_cellid_format = True

    rows = []
    for feat in data["features"]:
        cell_id = feat["properties"]["cell_id"]
        coords = np.array(feat["geometry"]["coordinates"][0])  # exterior ring
        cx, cy = coords.mean(axis=0)
        if use_cellid_format:
            idx = f"cellid_{int(cell_id):09d}-1"
        else:
            idx = str(cell_id)
        rows.append({"cell_id": idx, "x_centroid": cx, "y_centroid": cy})

    df = pd.DataFrame(rows).set_index("cell_id")
    return df


def load_visium_hd_cell_sample(
    spaceranger_dir: str | Path,
    sample_id: str,
    *,
    load_images: bool = False,
) -> ad.AnnData:
    """Load one Visium HD sample from SpaceRanger 4 cell segmentation output.

    Reads the cell-level (not bin-level) data produced by SpaceRanger 4.
    Searches for the segmentation output under both ``outs/segmented_outputs/``
    (SR4 default) and ``outs/cell_segmentation/`` (legacy). The expression
    matrix is ``filtered_feature_cell_matrix.h5`` (SR4) or
    ``filtered_feature_bc_matrix.h5`` (legacy). Cell spatial coordinates are
    extracted from ``cell_segmentations.geojson`` polygon centroids.

    Parameters
    ----------
    spaceranger_dir
        Path to Space Ranger output directory (containing ``outs/``).
    sample_id
        Sample identifier to store in obs['sample'].
    load_images
        Whether to load images.

    Returns
    -------
    AnnData with cell-level spatial coordinates and sample metadata.
    """
    import scanpy as sc

    spaceranger_dir = Path(spaceranger_dir)

    # Locate segmented output directory (SR4: segmented_outputs, legacy: cell_segmentation)
    cell_seg_dir = None
    for subdir in ["segmented_outputs", "cell_segmentation"]:
        for parent in [spaceranger_dir / "outs", spaceranger_dir]:
            candidate = parent / subdir
            if candidate.exists():
                cell_seg_dir = candidate
                break
        if cell_seg_dir is not None:
            break

    if cell_seg_dir is None:
        raise FileNotFoundError(
            f"Cell segmentation directory not found under {spaceranger_dir}. "
            f"Tried: outs/segmented_outputs, outs/cell_segmentation, "
            f"segmented_outputs, cell_segmentation"
        )

    # Load the filtered expression matrix (SR4: cell_matrix, legacy: bc_matrix)
    h5_candidates = [
        cell_seg_dir / "filtered_feature_cell_matrix.h5",
        cell_seg_dir / "filtered_feature_bc_matrix.h5",
    ]
    mtx_candidates = [
        cell_seg_dir / "filtered_feature_cell_matrix",
        cell_seg_dir / "filtered_feature_bc_matrix",
    ]

    h5_path = next((p for p in h5_candidates if p.exists()), None)
    mtx_dir = next((p for p in mtx_candidates if p.exists()), None)

    if h5_path is not None:
        adata = sc.read_10x_h5(str(h5_path))
    elif mtx_dir is not None:
        adata = sc.read_10x_mtx(str(mtx_dir))
    else:
        raise FileNotFoundError(
            f"No filtered matrix found in {cell_seg_dir}. "
            f"Expected filtered_feature_cell_matrix.h5 or filtered_feature_bc_matrix.h5"
        )

    # Load cell coordinates — try geojson first (SR4), then parquet/CSV fallbacks
    coords_loaded = False

    # Strategy 1: GeoJSON cell segmentation polygons -> centroids
    geojson_path = cell_seg_dir / "cell_segmentations.geojson"
    if geojson_path.exists():
        try:
            centroids = _extract_centroids_from_geojson(geojson_path, obs_names=adata.obs_names)
            common = adata.obs_names.intersection(centroids.index)
            if len(common) > 0:
                adata = adata[common].copy()
                adata.obsm["spatial"] = np.array(
                    centroids.loc[common, ["x_centroid", "y_centroid"]]
                )
                coords_loaded = True
                logger.info(
                    "Loaded %d cell centroids from geojson for %s",
                    len(common),
                    sample_id,
                )
        except Exception as exc:
            logger.warning(
                "Failed to extract centroids from geojson for %s: %s",
                sample_id,
                exc,
            )

    # Strategy 2: parquet / CSV coordinate files
    if not coords_loaded:
        for coords_name in ["cells.parquet", "cells.csv.gz", "cells.csv"]:
            coords_path = cell_seg_dir / coords_name
            if coords_path.exists():
                if coords_name.endswith(".parquet"):
                    coords = pd.read_parquet(coords_path)
                else:
                    coords = pd.read_csv(coords_path)

                x_col = next(
                    (
                        c
                        for c in ["x_centroid", "cell_centroid_x", "pxl_col_in_fullres"]
                        if c in coords.columns
                    ),
                    None,
                )
                y_col = next(
                    (
                        c
                        for c in ["y_centroid", "cell_centroid_y", "pxl_row_in_fullres"]
                        if c in coords.columns
                    ),
                    None,
                )

                if x_col and y_col:
                    bc_col = next((c for c in ["barcode", "cell_id"] if c in coords.columns), None)
                    if bc_col:
                        coords = coords.set_index(bc_col)
                        common = adata.obs_names.intersection(coords.index)
                        if len(common) > 0:
                            adata = adata[common].copy()
                            adata.obsm["spatial"] = np.array(coords.loc[common, [x_col, y_col]])
                            coords_loaded = True
                    elif len(coords) == adata.n_obs:
                        adata.obsm["spatial"] = np.array(coords[[x_col, y_col]])
                        coords_loaded = True
                break

    # Strategy 3: tissue_positions.parquet
    if not coords_loaded:
        pos_file = cell_seg_dir / "spatial" / "tissue_positions.parquet"
        if pos_file.exists():
            df_pos = pd.read_parquet(pos_file).set_index("barcode")
            common = adata.obs_names.intersection(df_pos.index)
            if len(common) > 0:
                adata = adata[common].copy()
                if "pxl_col_in_fullres" in df_pos.columns:
                    adata.obsm["spatial"] = np.array(
                        df_pos.loc[common, ["pxl_col_in_fullres", "pxl_row_in_fullres"]]
                    )
                    coords_loaded = True

    if not coords_loaded:
        logger.warning("Could not load spatial coordinates for %s cell segmentation", sample_id)

    adata.obs["sample"] = sample_id
    adata.obs["library_id"] = sample_id
    adata.obs["raw_data_dir"] = str(spaceranger_dir)
    adata.var_names_make_unique()

    logger.info(
        "Loaded Visium HD cell segmentation sample %s: %d cells x %d genes",
        sample_id,
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
    processed_dir: str | Path,
    sample_id: str,
    *,
    load_images: bool = False,
    panel_csv: str | Path | None = None,
    rgb_channels: tuple[str, str, str] = ("PanCK", "CD3", "DNA1"),
    image_downsample: int = 1,
) -> ad.AnnData:
    """Load one IMC sample from a steinbock/ElementoLab processed directory.

    Expects the standard ElementoLab IMC pipeline output layout::

        processed/{sample}/
            tiffs/
                {roi_id}_full.tiff     # (C, H, W) multi-channel TIFF stack
                {roi_id}_full.csv      # channel index -> MarkerName(IsotopeTag)
                {roi_id}_full_mask.tiff
                {roi_id}_Probabilities.tiff
            cells.h5ad                 # single-cell data (steinbock default)

    Parameters
    ----------
    processed_dir
        Path to the processed sample directory (e.g. ``processed/{sample}/``).
        Must contain ``cells.h5ad`` (steinbock default) or ``cells/cells.h5ad``.
    sample_id
        Sample/ROI identifier (e.g. ``05122023_Vivek_S21_5251_G3_Group1-01``).
        Used to locate ``tiffs/{sample_id}_full.tiff`` when ``load_images=True``.
    load_images
        If ``True``, read the ``{sample_id}_full.tiff`` multi-channel stack
        from ``processed_dir/tiffs/``, build an arcsinh-normalized full stack
        and RGB composite, and store both in ``adata.uns['spatial'][sample_id]``
        (squidpy/scanpy-compatible format).
    panel_csv
        Optional path to a panel metadata CSV (``channel_labels.csv``) with
        columns ``channel, Target, Metal_Tag, Atom, full, ilastik``. Provides
        additional name aliases and flags. When ``None``, channels are resolved
        from the per-ROI ``*_full.csv`` alone.
    rgb_channels
        Three marker names ``(R, G, B)`` for the default RGB composite image.
        Defaults to ``("PanCK", "CD3", "DNA1")``.
    image_downsample
        Spatial downsampling factor applied uniformly (1 = no downsampling).
        Use 2 or 4 to reduce memory for large IMC images.

    Returns
    -------
    AnnData with sample annotation and spatial coordinates. When
    ``load_images=True``, ``adata.uns['spatial'][sample_id]`` contains::

        images/
          hires   (H, W, 3) uint8    — percentile-clipped RGB composite
          full    (C, H, W) float32  — arcsinh(x/5) normalized full stack
        scalefactors/
          tissue_hires_scalef: 1.0/downsample
          spot_diameter_fullres: 1.0  (1 pixel ~ 1 µm in IMC)
        metadata/
          channels: list[str]        — ordered protein names
          channel_strings: list[str] — full MarkerName(IsotopeTag) strings
          rgb_channels: {R, G, B}    — resolved protein names used
          rgb_indices: {R, G, B}     — TIFF stack indices used
          pixel_size_um: 1.0
    """
    from .imc import IMCPanelMapper, build_imc_composite

    # Support remote URIs: use smart_read_h5ad when processed_dir looks like
    # a direct h5ad URI, otherwise treat as a local directory path.
    processed_dir_str = str(processed_dir)
    if processed_dir_str.endswith(".h5ad") or "://" in processed_dir_str:
        # Direct path to a cells.h5ad (local or remote URI)
        try:
            from sc_tools.storage import smart_read_h5ad

            adata = smart_read_h5ad(processed_dir_str)
        except ImportError:
            adata = ad.read_h5ad(processed_dir_str)
        processed_dir = Path(processed_dir_str).parent
    else:
        processed_dir = Path(processed_dir)
        # Locate the cells h5ad — try common steinbock output locations
        candidates = [
            processed_dir / "cells.h5ad",
            processed_dir / "cells" / "cells.h5ad",
        ]
        h5ad_path = next((p for p in candidates if p.exists()), None)
        if h5ad_path is None:
            raise FileNotFoundError(
                f"Could not find cells.h5ad in {processed_dir}. "
                f"Tried: {[str(p) for p in candidates]}"
            )
        adata = ad.read_h5ad(h5ad_path)
    adata.obs["sample"] = sample_id
    adata.obs["library_id"] = sample_id
    adata.obs["raw_data_dir"] = str(processed_dir)
    adata.var_names_make_unique()

    # Ensure spatial coordinates exist
    if "spatial" not in adata.obsm and "X_spatial" in adata.obsm:
        adata.obsm["spatial"] = adata.obsm["X_spatial"]

    # Optionally load TIFF stack images
    if load_images:
        tiff_dir = processed_dir / "tiffs"
        if not tiff_dir.exists():
            logger.warning(
                "load_imc_sample: load_images=True but tiffs/ not found at %s; skipping",
                tiff_dir,
            )
        else:
            # Locate the *_full.tiff and *_full.csv for this ROI/sample.
            # Naming convention: {sample_id}_full.tiff / {sample_id}_full.csv
            tiff_path = tiff_dir / f"{sample_id}_full.tiff"
            csv_path = tiff_dir / f"{sample_id}_full.csv"

            # Fallback: glob for any *_full.tiff if direct path not found
            if not tiff_path.exists():
                candidates = sorted(tiff_dir.glob("*_full.tiff"))
                if candidates:
                    tiff_path = candidates[0]
                    roi_stem = tiff_path.name[: -len("_full.tiff")]
                    csv_path = tiff_dir / f"{roi_stem}_full.csv"
                    logger.warning(
                        "load_imc_sample: exact TIFF %s not found; using %s",
                        f"{sample_id}_full.tiff",
                        tiff_path.name,
                    )

            if not tiff_path.exists():
                logger.warning(
                    "load_imc_sample: no *_full.tiff found in %s; skipping image load",
                    tiff_dir,
                )
            elif not csv_path.exists():
                logger.warning(
                    "load_imc_sample: channel CSV %s not found; skipping image load",
                    csv_path,
                )
            else:
                try:
                    mapper = IMCPanelMapper(full_csv=csv_path, panel_csv=panel_csv)
                    # Also register var_names so partial matching works without panel CSV
                    if mapper.n_channels() == 0:
                        mapper.set_from_var_names(list(adata.var_names))

                    spatial_dict = build_imc_composite(
                        tiff_path=tiff_path,
                        channel_csv=csv_path,
                        panel_mapper=mapper,
                        r=rgb_channels[0],
                        g=rgb_channels[1],
                        b=rgb_channels[2],
                        downsample=image_downsample,
                    )
                    if "spatial" not in adata.uns:
                        adata.uns["spatial"] = {}
                    adata.uns["spatial"][sample_id] = spatial_dict
                    logger.info(
                        "Loaded IMC images for %s: %d channels, composite R=%s G=%s B=%s",
                        sample_id,
                        len(spatial_dict["metadata"]["channels"]),
                        spatial_dict["metadata"]["rgb_channels"].get("R"),
                        spatial_dict["metadata"]["rgb_channels"].get("G"),
                        spatial_dict["metadata"]["rgb_channels"].get("B"),
                    )
                except Exception as exc:
                    logger.warning(
                        "load_imc_sample: image loading failed for %s: %s", sample_id, exc
                    )

    logger.info(
        "Loaded IMC sample %s: %d cells x %d markers",
        sample_id,
        adata.n_obs,
        adata.n_vars,
    )
    return adata


def load_he_image(
    he_path: str | Path,
    library_id: str,
    adata: ad.AnnData,
    *,
    downsample: int = 1,
    image_key: str = "hires",
) -> None:
    """Load an H&E TIFF and inject it into ``adata.uns['spatial'][library_id]``.

    Works for any modality (IMC, Xenium, Visium HD) that has an accompanying
    H&E TIFF image. Creates the spatial dict if absent; overwrites
    ``images[image_key]`` if the key already exists.

    Parameters
    ----------
    he_path
        Path to the H&E TIFF file (any format readable by ``tifffile``).
    library_id
        Key under ``adata.uns['spatial']`` where the image will be stored.
    adata
        AnnData to modify in place.
    downsample
        Spatial downsampling factor applied uniformly (1 = no downsampling).
    image_key
        Key within ``adata.uns['spatial'][library_id]['images']`` under which
        the image is stored (default ``"hires"``).
    """
    try:
        import tifffile
    except ImportError as e:
        raise ImportError("tifffile is required for H&E image loading: pip install tifffile") from e

    he_path = Path(he_path)
    if not he_path.exists():
        raise FileNotFoundError(f"H&E image not found: {he_path}")

    img = tifffile.imread(str(he_path))

    # Ensure (H, W, C) layout
    if img.ndim == 2:
        # Grayscale -> RGB
        img = np.stack([img, img, img], axis=-1)
    elif img.ndim == 3 and img.shape[0] in (3, 4):
        # (C, H, W) -> (H, W, C)
        img = np.moveaxis(img, 0, -1)

    # Take RGB only (drop alpha if present)
    img = img[..., :3]

    if downsample > 1:
        img = img[::downsample, ::downsample, :]

    # Ensure uint8
    if img.dtype != np.uint8:
        if img.max() <= 1.0:
            img = (img * 255).clip(0, 255).astype(np.uint8)
        else:
            img = img.clip(0, 255).astype(np.uint8)

    scalef = 1.0 / downsample

    if "spatial" not in adata.uns:
        adata.uns["spatial"] = {}
    if library_id not in adata.uns["spatial"]:
        adata.uns["spatial"][library_id] = {"images": {}, "scalefactors": {}}

    spatial = adata.uns["spatial"][library_id]
    if "images" not in spatial:
        spatial["images"] = {}
    if "scalefactors" not in spatial:
        spatial["scalefactors"] = {}

    spatial["images"][image_key] = img
    spatial["scalefactors"].setdefault("tissue_hires_scalef", scalef)
    spatial["scalefactors"].setdefault("tissue_lowres_scalef", scalef)
    spatial["scalefactors"].setdefault("spot_diameter_fullres", 1.0)

    logger.info(
        "Loaded H&E image for library %s: shape %s (key=%s)",
        library_id,
        img.shape,
        image_key,
    )


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
