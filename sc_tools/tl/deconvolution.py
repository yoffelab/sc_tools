"""
Cell-type deconvolution for spatial transcriptomics.

Provides a generic ``deconvolution()`` function that dispatches to pluggable
backends (cell2location, tangram, destvi).  Per-library processing with backed
AnnData loading keeps peak memory low.

Typical usage
-------------
>>> import sc_tools as st
>>> st.tl.deconvolution(
...     spatial_adata="results/adata.normalized.scored.p35.h5ad",
...     sc_adata="results/seurat_object.h5ad",
...     method="cell2location",
...     celltype_key="cell.type",
... )

Existing helper ``select_signature_genes`` is preserved unchanged.
"""

from __future__ import annotations

import logging
import pickle
from pathlib import Path
from typing import TYPE_CHECKING, Any, Protocol, runtime_checkable

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from ..data.io import get_cache_key, load_cached_signatures, save_cached_signatures
from ..memory.gpu import get_gpu_setting
from ..memory.profiling import (
    aggressive_cleanup,
    check_memory_threshold,
    log_memory,
)

if TYPE_CHECKING:
    from logging import Logger

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Backend protocol and registry
# ---------------------------------------------------------------------------


@runtime_checkable
class DeconvolutionBackend(Protocol):
    """Protocol every deconvolution backend must satisfy."""

    @staticmethod
    def run(
        sc_adata: ad.AnnData,
        spatial_adata_lib: ad.AnnData,
        shared_genes: list[str],
        celltype_key: str,
        *,
        use_gpu: bool = False,
        reference_profiles: pd.DataFrame | None = None,
        logger_instance: Logger | None = None,
        **kwargs: Any,
    ) -> np.ndarray | None:
        """Return an (n_spots, n_celltypes) proportion matrix or *None* on failure."""
        ...


_BACKENDS: dict[str, type[DeconvolutionBackend]] = {}


def register_backend(name: str, cls: type[DeconvolutionBackend]) -> None:
    """Register a backend class under *name*."""
    _BACKENDS[name] = cls


def get_backend(name: str) -> type[DeconvolutionBackend]:
    """Retrieve a registered backend by *name*."""
    if name not in _BACKENDS:
        available = ", ".join(sorted(_BACKENDS)) or "(none)"
        raise ValueError(f"Unknown deconvolution method {name!r}. Available: {available}")
    return _BACKENDS[name]


# ---------------------------------------------------------------------------
# Reference profile extraction (memory optimisation)
# ---------------------------------------------------------------------------


def extract_reference_profiles(
    sc_adata: ad.AnnData,
    celltype_key: str,
    genes: list[str] | None = None,
    qc_labels: list[str] | None = None,
    cache_path: str | Path | None = None,
    logger_instance: Logger | None = None,
) -> pd.DataFrame:
    """Compute mean expression per cell type from scRNA-seq reference.

    The resulting DataFrame (genes x cell_types) is ~100x smaller than the
    full reference and can be passed directly to Cell2location via its
    ``cell_state_df`` parameter, skipping regression model training.

    Parameters
    ----------
    sc_adata
        Single-cell reference AnnData.
    celltype_key
        Column in ``sc_adata.obs`` with cell-type labels.
    genes
        Subset of genes to include.  *None* keeps all.
    qc_labels
        Cell-type labels to exclude (e.g. ``["Doublets", "QC_Filtered"]``).
    cache_path
        If given, cache the result as a pickle file.
    logger_instance
        Optional logger.

    Returns
    -------
    pandas.DataFrame
        Genes (rows) x cell types (columns) mean expression matrix.
    """
    log = logger_instance or logger

    # Try loading from cache first
    if cache_path is not None:
        cache_path = Path(cache_path)
        if cache_path.exists():
            try:
                with open(cache_path, "rb") as fh:
                    cached = pickle.load(fh)
                if isinstance(cached, pd.DataFrame):
                    log.info(
                        f"Loaded reference profiles from cache: {cache_path} "
                        f"({cached.shape[0]} genes x {cached.shape[1]} cell types)"
                    )
                    return cached
            except Exception as exc:
                log.warning(f"Failed to load cached reference profiles: {exc}")

    adata = sc_adata
    if qc_labels:
        mask = ~adata.obs[celltype_key].isin(qc_labels)
        adata = adata[mask]
        log.info(f"Excluded {(~mask).sum()} cells matching QC labels")

    if genes is not None:
        shared = [g for g in genes if g in adata.var_names]
        adata = adata[:, shared]
        log.info(f"Subset to {len(shared)} genes")

    expr_adata = adata
    expr_var_names = adata.var_names

    log.info("Computing mean expression per cell type...")
    celltypes = adata.obs[celltype_key]
    unique_cts = sorted(celltypes.unique())
    profiles: dict[str, np.ndarray] = {}
    for ct in unique_cts:
        ct_mask = celltypes == ct
        X_ct = expr_adata[ct_mask].X
        if hasattr(X_ct, "toarray"):
            X_ct = X_ct.toarray()
        profiles[ct] = np.asarray(X_ct, dtype=np.float64).mean(axis=0)

    df = pd.DataFrame(profiles, index=expr_var_names)

    # Ensure strictly positive values (required for Cell2location GammaPoisson)
    min_val = float(df.values.min())
    if min_val < 0:
        # SCTransform Pearson residuals can be negative. Shift all values so minimum
        # is positive, preserving relative differences between cell types.
        shift = abs(min_val) + 1.0
        log.info(
            f"Reference profiles have negative values (min={min_val:.2f}); "
            f"shifting by +{shift:.2f} for Cell2location compatibility"
        )
        df = df + shift

    # Add small epsilon to avoid exact zeros (GammaPoisson rate must be > 0)
    df = df.clip(lower=1e-10)

    log.info(f"Reference profiles: {df.shape[0]} genes x {df.shape[1]} cell types")

    if cache_path is not None:
        try:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            with open(cache_path, "wb") as fh:
                pickle.dump(df, fh)
            log.info(f"Cached reference profiles to {cache_path}")
        except Exception as exc:
            log.warning(f"Failed to cache reference profiles: {exc}")

    return df


# ---------------------------------------------------------------------------
# Backend: Tangram
# ---------------------------------------------------------------------------


class TangramBackend:
    """Tangram (OT-based) deconvolution backend."""

    @staticmethod
    def run(
        sc_adata: ad.AnnData,
        spatial_adata_lib: ad.AnnData,
        shared_genes: list[str],
        celltype_key: str,
        *,
        use_gpu: bool = False,
        reference_profiles: pd.DataFrame | None = None,
        logger_instance: Logger | None = None,
        **kwargs: Any,
    ) -> np.ndarray | None:
        log = logger_instance or logger
        try:
            import tangram as tg
        except ImportError:
            log.warning("tangram-sc not installed, skipping Tangram backend")
            return None

        num_epochs = kwargs.get("num_epochs", 500)
        log.info(f"Tangram: mapping with {num_epochs} epochs on {len(shared_genes)} genes")

        try:
            sc_copy = sc_adata[:, shared_genes].copy()
            sp_copy = spatial_adata_lib[:, shared_genes].copy()
            log_memory("Tangram: after gene subset", sp_copy, logger_instance=log)

            if not check_memory_threshold(
                threshold_mb=50000, threshold_percent=92.0, logger_instance=log
            ):
                log.warning("Memory too high for Tangram, aborting")
                del sc_copy, sp_copy
                aggressive_cleanup()
                return None

            tg.pp_adatas(sc_copy, sp_copy, genes=shared_genes)

            # Safety: check for NaN/Inf in inputs before Tangram mapping
            for _label, _ad in [("sc", sc_copy), ("spatial", sp_copy)]:
                _Xc = _ad.X
                if hasattr(_Xc, "toarray"):
                    _Xc = _Xc.toarray()
                _n_nan = int(np.isnan(_Xc).sum())
                if _n_nan > 0:
                    log.warning(
                        f"Tangram: {_label} has {_n_nan} NaN values in X — replacing with 0"
                    )
                    _Xc = np.nan_to_num(_Xc, nan=0.0, posinf=0.0, neginf=0.0)
                    import scipy.sparse as _spt

                    if _spt.issparse(_ad.X):
                        _ad.X = _spt.csr_matrix(_Xc)
                    else:
                        _ad.X = _Xc

            ad_map = tg.map_cells_to_space(
                sc_copy,
                sp_copy,
                mode="clusters",
                cluster_label=celltype_key,
                num_epochs=num_epochs,
            )

            # Check if mapping produced NaN (can happen with problematic data)
            _map_X = ad_map.X
            if hasattr(_map_X, "toarray"):
                _map_X = _map_X.toarray()
            if np.isnan(_map_X).all():
                log.error(
                    "Tangram: mapping matrix is entirely NaN — "
                    "likely due to data quality issues in reference or spatial data"
                )
                del sc_copy, sp_copy, ad_map
                aggressive_cleanup()
                return None

            del sc_copy, sp_copy
            aggressive_cleanup()

            # Extract proportions using project_cell_annotations
            tg.project_cell_annotations(ad_map, spatial_adata_lib, annotation=celltype_key)

            proportions = _extract_tangram_proportions(
                ad_map, spatial_adata_lib, sc_adata, celltype_key, log
            )
            del ad_map
            aggressive_cleanup()
            return proportions

        except Exception as exc:
            log.error(f"Tangram failed: {exc}")
            aggressive_cleanup()
            return None


def _extract_tangram_proportions(
    ad_map: ad.AnnData,
    spatial_adata_lib: ad.AnnData,
    sc_adata: ad.AnnData,
    celltype_key: str,
    log: Logger,
) -> np.ndarray | None:
    """Extract cell-type proportions from Tangram output.

    Handles both sparse and dense matrices, auto-detects orientation
    (spots x celltypes vs celltypes x spots).  Ported from
    robin ``ct_proportions_from_tangram``.
    """
    import scipy.sparse as sp

    unique_cts = sorted(sc_adata.obs[celltype_key].unique())

    # First try obsm keys set by project_cell_annotations
    for key in spatial_adata_lib.obsm:
        if "proportion" in key.lower() or celltype_key.lower() in key.lower():
            prop_data = spatial_adata_lib.obsm[key]
            if hasattr(prop_data, "values"):
                prop_data = prop_data.values
            if prop_data.shape[0] == spatial_adata_lib.n_obs:
                log.info(f"Tangram: extracted proportions from obsm[{key!r}]")
                return np.asarray(prop_data, dtype=np.float32)

    # Try obs columns matching cell type names
    ct_cols = [ct for ct in unique_cts if ct in spatial_adata_lib.obs.columns]
    if len(ct_cols) >= len(unique_cts) * 0.8:
        log.info(f"Tangram: extracting proportions from {len(ct_cols)} obs columns")
        parts = []
        for ct in unique_cts:
            if ct in spatial_adata_lib.obs.columns:
                parts.append(spatial_adata_lib.obs[ct].values)
            else:
                parts.append(np.zeros(spatial_adata_lib.n_obs))
        proportions = np.column_stack(parts).astype(np.float32)
        row_sums = proportions.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        return proportions / row_sums

    # Fallback: normalise the mapping matrix itself (robin-style)
    X = ad_map.X
    if sp.issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float32)

    # Auto-detect orientation: more cols than rows -> spots are columns
    if X.shape[1] > X.shape[0]:
        s = X.sum(axis=0)
        s[s == 0] = 1.0
        proportions = (X / s).T
    else:
        s = X.sum(axis=1, keepdims=True)
        s[s == 0] = 1.0
        proportions = X / s

    log.info(f"Tangram: extracted proportions from mapping matrix ({proportions.shape})")
    return proportions


# ---------------------------------------------------------------------------
# Backend: Cell2location
# ---------------------------------------------------------------------------


class Cell2locationBackend:
    """Cell2location deconvolution backend (GPU-recommended)."""

    @staticmethod
    def run(
        sc_adata: ad.AnnData,
        spatial_adata_lib: ad.AnnData,
        shared_genes: list[str],
        celltype_key: str,
        *,
        use_gpu: bool = False,
        reference_profiles: pd.DataFrame | None = None,
        logger_instance: Logger | None = None,
        **kwargs: Any,
    ) -> np.ndarray | None:
        log = logger_instance or logger
        try:
            import cell2location as c2l
            from cell2location.models import Cell2location, RegressionModel
        except ImportError:
            log.warning("cell2location not installed, skipping")
            return None

        # CPU-aware defaults: 25000 epochs on GPU, 1000 on CPU
        default_epochs = 25000 if use_gpu else 1000
        max_epochs = kwargs.get("max_epochs", default_epochs)
        num_samples = kwargs.get("num_samples", 1000)
        batch_size = kwargs.get("batch_size", 2500)
        if not use_gpu:
            log.info(f"Cell2location: running on CPU (max_epochs={max_epochs})")

        try:
            log_memory("Cell2location: start", logger_instance=log)

            mem_threshold_mb = kwargs.get("memory_threshold_mb", 50000)
            mem_threshold_pct = kwargs.get("memory_threshold_pct", 92.0)
            if not check_memory_threshold(
                threshold_mb=mem_threshold_mb,
                threshold_percent=mem_threshold_pct,
                logger_instance=log,
            ):
                log.warning("Memory too high for Cell2location, aborting")
                return None

            spatial_lib_sig = spatial_adata_lib[:, shared_genes].copy()

            if spatial_lib_sig.raw is None:
                spatial_lib_sig.raw = spatial_lib_sig.copy()

            # --- Memory-optimised path: use pre-computed reference profiles ---
            if reference_profiles is not None:
                log.info(
                    "Cell2location: using pre-computed reference profiles (skipping regression)"
                )
                # Align genes
                common_genes = [g for g in shared_genes if g in reference_profiles.index]
                cell_state_df = reference_profiles.loc[common_genes]
                spatial_lib_sig = spatial_lib_sig[:, common_genes].copy()
                if spatial_lib_sig.raw is None:
                    spatial_lib_sig.raw = spatial_lib_sig.copy()

                c2l.models.Cell2location.setup_anndata(spatial_lib_sig, batch_key=None)
                mod_spatial = Cell2location(
                    spatial_lib_sig,
                    cell_state_df=cell_state_df,
                )
            else:
                # --- Full path: train regression model on reference ---
                log.info("Cell2location: training regression model on reference")
                sc_copy = sc_adata[:, shared_genes].copy()
                if sc_copy.raw is None:
                    sc_copy.raw = sc_copy.copy()

                c2l.models.RegressionModel.setup_anndata(
                    sc_copy, labels_key=celltype_key, batch_key=None
                )
                mod = RegressionModel(sc_copy)

                # Adaptive epoch reduction under memory pressure
                actual_epochs = max_epochs
                if not check_memory_threshold(
                    threshold_mb=50000, threshold_percent=90.0, logger_instance=log
                ):
                    actual_epochs = min(max_epochs, 15000)
                    num_samples = min(num_samples, 500)
                    batch_size = min(batch_size, 1500)
                    log.info(
                        f"Cell2location: reducing parameters (epochs={actual_epochs}) due to memory"
                    )

                mod.train(max_epochs=actual_epochs)
                sc_copy = mod.export_posterior(
                    sc_copy,
                    sample_kwargs={"num_samples": num_samples, "batch_size": batch_size},
                )

                c2l.models.Cell2location.setup_anndata(spatial_lib_sig, batch_key=None)
                mod_spatial = Cell2location(spatial_lib_sig, sc_copy)
                del sc_copy, mod
                aggressive_cleanup()

            # Adaptive epoch reduction
            actual_epochs = max_epochs
            if not check_memory_threshold(
                threshold_mb=50000, threshold_percent=90.0, logger_instance=log
            ):
                actual_epochs = min(max_epochs, 15000)
                log.info(f"Cell2location: reducing spatial epochs to {actual_epochs}")

            mod_spatial.train(max_epochs=actual_epochs)
            spatial_lib_sig = mod_spatial.export_posterior(
                spatial_lib_sig,
                sample_kwargs={"num_samples": num_samples, "batch_size": batch_size},
            )

            # Extract proportions
            proportions = _extract_c2l_proportions(spatial_lib_sig, log)

            del mod_spatial, spatial_lib_sig
            aggressive_cleanup()
            log_memory("Cell2location: done", logger_instance=log)
            return proportions

        except Exception as exc:
            log.error(f"Cell2location failed: {exc}")
            aggressive_cleanup()
            return None


def _extract_c2l_proportions(spatial_lib_sig: ad.AnnData, log: Logger) -> np.ndarray | None:
    """Normalise Cell2location abundance to row-sum-one proportions."""
    for key in ("q05_cell_abundance_w_sf", "means_cell_abundance_w_sf"):
        if key in spatial_lib_sig.obsm:
            vals = spatial_lib_sig.obsm[key]
            if hasattr(vals, "values"):
                vals = vals.values
            proportions = np.asarray(vals, dtype=np.float64).copy()
            row_sums = proportions.sum(axis=1, keepdims=True)
            row_sums[row_sums == 0] = 1.0
            proportions = (proportions / row_sums).astype(np.float32)
            log.info(f"Cell2location: proportions from obsm[{key!r}] ({proportions.shape})")
            return proportions

    log.warning(
        f"Cell2location: no proportion key found in obsm. "
        f"Available: {list(spatial_lib_sig.obsm.keys())}"
    )
    return None


# ---------------------------------------------------------------------------
# Backend: DestVI
# ---------------------------------------------------------------------------


class DestVIBackend:
    """DestVI (scvi-tools) deconvolution backend."""

    @staticmethod
    def run(
        sc_adata: ad.AnnData,
        spatial_adata_lib: ad.AnnData,
        shared_genes: list[str],
        celltype_key: str,
        *,
        use_gpu: bool = False,
        reference_profiles: pd.DataFrame | None = None,
        logger_instance: Logger | None = None,
        **kwargs: Any,
    ) -> np.ndarray | None:
        log = logger_instance or logger
        try:
            from scvi.external import DestVI
        except ImportError:
            log.warning("scvi-tools / DestVI not installed, skipping")
            return None

        max_epochs = kwargs.get("max_epochs", 25000)

        try:
            log_memory("DestVI: start", logger_instance=log)

            if not check_memory_threshold(
                threshold_mb=50000, threshold_percent=92.0, logger_instance=log
            ):
                log.warning("Memory too high for DestVI, aborting")
                return None

            sc_copy = sc_adata[:, shared_genes].copy()
            spatial_lib_sig = spatial_adata_lib[:, shared_genes].copy()

            if sc_copy.raw is None:
                sc_copy.raw = sc_copy.copy()
            if spatial_lib_sig.raw is None:
                spatial_lib_sig.raw = spatial_lib_sig.copy()

            # Adaptive epochs
            actual_epochs = max_epochs
            if not check_memory_threshold(
                threshold_mb=50000, threshold_percent=90.0, logger_instance=log
            ):
                actual_epochs = min(max_epochs, 10000)
                log.info(f"DestVI: reducing epochs to {actual_epochs}")

            DestVI.setup_anndata(sc_copy, labels_key=celltype_key)
            vae = DestVI(sc_copy)
            vae.train(max_epochs=actual_epochs)

            del sc_copy
            aggressive_cleanup()

            DestVI.setup_anndata(spatial_lib_sig, labels_key=None)
            vae_st = DestVI.load_query_data(spatial_lib_sig, vae)
            vae_st.train(max_epochs=actual_epochs)

            proportions_df = vae_st.get_proportions(spatial_lib_sig)
            proportions = np.asarray(proportions_df.values, dtype=np.float32)
            log.info(f"DestVI: proportions ({proportions.shape})")

            del vae, vae_st, spatial_lib_sig
            aggressive_cleanup()
            return proportions

        except Exception as exc:
            log.error(f"DestVI failed: {exc}")
            aggressive_cleanup()
            return None


# ---------------------------------------------------------------------------
# Register built-in backends
# ---------------------------------------------------------------------------

register_backend("tangram", TangramBackend)
register_backend("cell2location", Cell2locationBackend)
register_backend("destvi", DestVIBackend)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def deconvolution(
    spatial_adata: ad.AnnData | str | Path,
    sc_adata: ad.AnnData | str | Path | None = None,
    *,
    method: str = "cell2location",
    celltype_key: str = "celltype",
    spatial_batch_key: str = "library_id",
    sc_batch_key: str | None = None,
    reference_profiles: pd.DataFrame | Path | None = None,
    n_signature_genes: int = 2000,
    use_gpu: bool | None = None,
    qc_labels: list[str] | None = None,
    method_kwargs: dict | None = None,
    cache_dir: str | Path | None = None,
    output_dir: str | Path | None = None,
    output_file: str | Path | None = None,
    logger_instance: Logger | None = None,
) -> ad.AnnData:
    """Run cell-type deconvolution on spatial transcriptomics data.

    Processes each library (``spatial_batch_key``) independently using backed
    AnnData loading to keep peak memory low.

    Parameters
    ----------
    spatial_adata
        Spatial AnnData or path to h5ad file.
    sc_adata
        Single-cell reference AnnData or path.  Required unless
        *reference_profiles* is provided (Cell2location only).
    method
        Backend name: ``"cell2location"`` (default), ``"tangram"``, ``"destvi"``.
    celltype_key
        Column in ``sc_adata.obs`` with cell-type labels.
    spatial_batch_key
        Column in spatial ``obs`` for per-library batching.
    sc_batch_key
        Batch key in ``sc_adata.obs`` (used for HVG selection).
    reference_profiles
        Pre-computed reference profiles (DataFrame or path to pickle).
        Memory optimisation for Cell2location -- skips regression training.
    n_signature_genes
        Number of signature genes for deconvolution.
    use_gpu
        GPU setting.  ``None`` = auto-detect.
    qc_labels
        Cell-type labels to exclude from the reference.
    method_kwargs
        Extra keyword arguments forwarded to the backend ``run()`` method.
    cache_dir
        Directory for caching signature genes and reference profiles.
    output_dir
        Directory for per-library intermediate results.
    output_file
        Path to save the final AnnData with proportions.
    logger_instance
        Optional logger.

    Returns
    -------
    AnnData
        The spatial AnnData with ``obsm['cell_type_proportions']`` (n_spots x
        n_celltypes) and ``obs['{method}_argmax']`` (dominant cell type per
        spot).
    """
    log = logger_instance or logger
    method_kwargs = method_kwargs or {}
    backend_cls = get_backend(method)

    # --- Resolve GPU ---
    gpu = get_gpu_setting(use_gpu)
    log.info(f"Deconvolution: method={method}, use_gpu={gpu}")

    # --- Resolve spatial data ---
    spatial_path: str | None = None
    if isinstance(spatial_adata, (str, Path)):
        spatial_path = str(spatial_adata)
        log.info(f"Loading spatial data from {spatial_path}")
        spatial_adata_full = sc.read_h5ad(spatial_path)
    else:
        spatial_adata_full = spatial_adata

    log_memory("After loading spatial data", spatial_adata_full, logger_instance=log)

    # --- Resolve scRNA-seq reference ---
    sc_adata_obj: ad.AnnData | None = None
    sc_data_file: str | None = None

    if isinstance(sc_adata, (str, Path)):
        sc_data_file = str(sc_adata)
        log.info(f"Loading scRNA-seq reference from {sc_data_file}")
        sc_adata_obj = sc.read_h5ad(sc_data_file)
        sc_adata_obj.var_names_make_unique()
        log_memory("After loading scRNA-seq", sc_adata_obj, logger_instance=log)
    elif sc_adata is not None:
        sc_adata_obj = sc_adata
    elif reference_profiles is None:
        raise ValueError("Either sc_adata or reference_profiles must be provided")

    # --- Preprocess reference ---
    unique_cts: list[str] = []
    if sc_adata_obj is not None:
        # Remove cells with zero total counts (they produce NaN after log1p)
        import scipy.sparse as _sp

        _X = sc_adata_obj.X
        if _sp.issparse(_X):
            total_counts = np.asarray(_X.sum(axis=1)).ravel()
        else:
            total_counts = np.asarray(_X.sum(axis=1)).ravel()
        zero_mask = total_counts == 0
        if zero_mask.any():
            n_zero = int(zero_mask.sum())
            log.warning(
                f"Removing {n_zero} cells with zero total counts from reference "
                f"(would produce NaN after normalisation)"
            )
            sc_adata_obj = sc_adata_obj[~zero_mask].copy()

        # Normalise if raw counts (skip if data has negative values — already scaled)
        if _sp.issparse(sc_adata_obj.X):
            _X_check = sc_adata_obj.X.toarray()
        else:
            _X_check = np.asarray(sc_adata_obj.X)
        _min_val = float(np.nanmin(_X_check))
        _max_val = float(np.nanmax(_X_check))
        del _X_check

        if _min_val < 0:
            log.info(
                f"Reference X has negative values (min={_min_val:.2f}) — "
                "already normalised/scaled, skipping normalisation"
            )
        elif _max_val > 100:
            log.info("Normalising scRNA-seq reference (raw counts detected)")
            sc.pp.normalize_total(sc_adata_obj, target_sum=1e4)
            sc.pp.log1p(sc_adata_obj)

        # Safety net: replace any remaining NaN/Inf in X with 0
        _X_post = sc_adata_obj.X
        if _sp.issparse(_X_post):
            _X_dense = _X_post.toarray()
        else:
            _X_dense = np.asarray(_X_post)
        n_nan = int(np.isnan(_X_dense).sum())
        n_inf = int(np.isinf(_X_dense).sum())
        if n_nan > 0 or n_inf > 0:
            log.warning(f"Replacing {n_nan} NaN and {n_inf} Inf values in reference X with 0")
            _X_dense = np.nan_to_num(_X_dense, nan=0.0, posinf=0.0, neginf=0.0)
            if _sp.issparse(sc_adata_obj.X):
                import scipy.sparse as sp2

                sc_adata_obj.X = sp2.csr_matrix(_X_dense)
            else:
                sc_adata_obj.X = _X_dense
        del _X_dense

        # Filter QC labels from cell-type list (keep cells for marker computation)
        all_cts = sorted(sc_adata_obj.obs[celltype_key].unique())
        if qc_labels:
            unique_cts = [ct for ct in all_cts if ct not in qc_labels]
        else:
            unique_cts = all_cts
        log.info(f"Cell types for deconvolution: {len(unique_cts)}")

    # --- Resolve reference profiles ---
    ref_profiles_df: pd.DataFrame | None = None
    if isinstance(reference_profiles, (str, Path)):
        with open(reference_profiles, "rb") as fh:
            ref_profiles_df = pickle.load(fh)
        log.info(f"Loaded reference profiles from {reference_profiles}")
    elif isinstance(reference_profiles, pd.DataFrame):
        ref_profiles_df = reference_profiles

    # --- Auto-compute reference profiles for cell2location ---
    # Cell2location's full path trains an NB regression model on all reference cells,
    # which needs raw integer counts and is very memory/time-intensive on CPU.
    # We auto-switch to the cell_state_df shortcut (mean expression profiles) when:
    #   1. The reference lacks raw integer counts (already normalised/scaled), OR
    #   2. Running on CPU (training regression on 100K+ cells on CPU is impractical)
    if method == "cell2location" and ref_profiles_df is None and sc_adata_obj is not None:
        _has_raw_counts = False
        # Check if .raw has integer counts
        if sc_adata_obj.raw is not None:
            _raw_sample = sc_adata_obj.raw.X[:100]
            if hasattr(_raw_sample, "toarray"):
                _raw_sample = _raw_sample.toarray()
            _raw_sample = np.asarray(_raw_sample)
            if np.allclose(_raw_sample, np.round(_raw_sample)) and _raw_sample.min() >= 0:
                _has_raw_counts = True
        # Check if .X has integer counts
        if not _has_raw_counts:
            _x_sample = sc_adata_obj.X[:100]
            if hasattr(_x_sample, "toarray"):
                _x_sample = _x_sample.toarray()
            _x_sample = np.asarray(_x_sample)
            if (
                np.allclose(_x_sample, np.round(_x_sample))
                and _x_sample.min() >= 0
                and _x_sample.max() > 100
            ):
                _has_raw_counts = True

        _use_profiles = False
        if not _has_raw_counts:
            log.info(
                "Cell2location: reference lacks raw integer counts — "
                "using reference profiles (cell_state_df) shortcut"
            )
            _use_profiles = True
        elif not gpu:
            log.info(
                "Cell2location: CPU mode — using reference profiles (cell_state_df) "
                "shortcut to avoid expensive NB regression training"
            )
            _use_profiles = True

        if _use_profiles:
            _cache_path = Path(cache_dir) / "reference_profiles.pkl" if cache_dir else None
            ref_profiles_df = extract_reference_profiles(
                sc_adata_obj,
                celltype_key,
                qc_labels=qc_labels,
                cache_path=_cache_path,
                logger_instance=log,
            )

    # --- Select signature genes ---
    candidate_genes: list[str] = []
    if sc_adata_obj is not None:
        candidate_genes = select_signature_genes(
            sc_adata_obj,
            celltype_key,
            sc_batch_key or "",
            n_signature_genes,
            skip_hvg=(sc_batch_key is None),
            cache_dir=str(cache_dir) if cache_dir else None,
            sc_data_file=sc_data_file,
            logger_instance=log,
        )
    elif ref_profiles_df is not None:
        candidate_genes = list(ref_profiles_df.index)

    # --- Free scRNA-seq reference if not needed for per-library processing ---
    # When using reference_profiles (cell2location cell_state_df path), the backend
    # only needs the profiles DataFrame, not the full scRNA-seq object.
    # For Tangram, the full reference IS needed (OT mapping uses all cells).
    if ref_profiles_df is not None and method == "cell2location" and sc_adata_obj is not None:
        log.info("Freeing scRNA-seq reference (reference_profiles available for Cell2location)")
        del sc_adata_obj
        sc_adata_obj = None
        aggressive_cleanup()
        log_memory("After freeing scRNA-seq reference", logger_instance=log)

    # --- Determine libraries ---
    if spatial_batch_key not in spatial_adata_full.obs.columns:
        # Fallback: try sample_id, or treat whole dataset as one batch
        if "sample_id" in spatial_adata_full.obs.columns:
            spatial_batch_key = "sample_id"
            log.info("Using 'sample_id' as spatial_batch_key")
        else:
            log.info("No batch key found; processing entire dataset as a single batch")
            spatial_adata_full.obs["_single_batch"] = "all"
            spatial_batch_key = "_single_batch"

    library_ids = sorted(spatial_adata_full.obs[spatial_batch_key].unique())
    log.info(f"Libraries to process: {len(library_ids)}: {library_ids}")

    # --- Per-library processing ---
    all_proportions: list[np.ndarray] = []
    all_spot_indices: list[str] = []
    n_celltypes: int | None = None

    if output_dir is not None:
        Path(output_dir).mkdir(parents=True, exist_ok=True)

    for lib_id in library_ids:
        log.info(f"{'=' * 60}")
        log.info(f"Processing library: {lib_id}")
        log_memory(f"Start {lib_id}", logger_instance=log)

        try:
            # Load library subset via backed file if path available
            if spatial_path is not None:
                backed = ad.read_h5ad(spatial_path, backed="r")
                lib_mask = backed.obs[spatial_batch_key] == lib_id
                lib_obs_names = backed.obs[lib_mask].index.tolist()
                spatial_adata_lib = backed[lib_mask].to_memory()
                del backed
            else:
                lib_mask = spatial_adata_full.obs[spatial_batch_key] == lib_id
                lib_obs_names = spatial_adata_full.obs[lib_mask].index.tolist()
                spatial_adata_lib = spatial_adata_full[lib_mask].copy()

            aggressive_cleanup()

            if len(lib_obs_names) == 0:
                log.warning(f"No spots for library {lib_id}, skipping")
                continue

            log.info(f"Library {lib_id}: {spatial_adata_lib.shape}")

            # Use raw counts if .raw exists
            if spatial_adata_lib.raw is not None:
                log.info("Using raw counts from spatial data")
                raw_ad = spatial_adata_lib.raw.to_adata()
                spatial_adata_lib.X = raw_ad[:, spatial_adata_lib.var_names].X
                del raw_ad
                aggressive_cleanup()

            # Find shared genes
            shared_genes = list(set(candidate_genes).intersection(spatial_adata_lib.var_names))
            if sc_adata_obj is not None:
                shared_genes = [g for g in shared_genes if g in sc_adata_obj.var_names]
            shared_genes = shared_genes[:n_signature_genes]
            log.info(f"Shared signature genes: {len(shared_genes)}")

            if len(shared_genes) < 100:
                log.warning(f"Too few shared genes ({len(shared_genes)}), skipping {lib_id}")
                continue

            # Reduce genes if memory is tight
            if (
                not check_memory_threshold(
                    threshold_mb=50000, threshold_percent=90.0, logger_instance=log
                )
                and len(shared_genes) > 1000
            ):
                shared_genes = shared_genes[:1000]
                log.info(f"Reduced to {len(shared_genes)} genes due to memory pressure")
                aggressive_cleanup()

            # Run backend
            proportions = backend_cls.run(
                sc_adata_obj if sc_adata_obj is not None else ad.AnnData(),
                spatial_adata_lib,
                shared_genes,
                celltype_key,
                use_gpu=gpu,
                reference_profiles=ref_profiles_df,
                logger_instance=log,
                **method_kwargs,
            )

            if proportions is not None and proportions.shape[0] == len(lib_obs_names):
                all_proportions.append(proportions)
                all_spot_indices.extend(lib_obs_names)
                if n_celltypes is None:
                    n_celltypes = proportions.shape[1]
                log.info(f"Library {lib_id}: proportions {proportions.shape}")
            elif proportions is not None:
                log.warning(
                    f"Shape mismatch for {lib_id}: "
                    f"proportions {proportions.shape[0]} vs spots {len(lib_obs_names)}"
                )
            else:
                log.warning(f"No proportions for library {lib_id}")

            del spatial_adata_lib
            aggressive_cleanup()
            log_memory(f"Done {lib_id}", logger_instance=log)

        except Exception as exc:
            log.error(f"Error processing library {lib_id}: {exc}")
            aggressive_cleanup()
            continue

    # --- Assemble proportions ---
    if len(all_proportions) == 0:
        log.warning("No proportions extracted from any library")
        return spatial_adata_full

    proportions_matrix = np.vstack(all_proportions)
    log.info(f"Total proportions: {proportions_matrix.shape}")

    # Align to spatial_adata_full order
    index_map = {idx: i for i, idx in enumerate(all_spot_indices)}
    ordered = np.zeros((spatial_adata_full.n_obs, proportions_matrix.shape[1]), dtype=np.float32)
    for j, obs_name in enumerate(spatial_adata_full.obs_names):
        if obs_name in index_map:
            ordered[j] = proportions_matrix[index_map[obs_name]]

    # Determine column names for the proportions DataFrame
    if unique_cts and len(unique_cts) == ordered.shape[1]:
        ct_names = unique_cts
    else:
        ct_names = [f"celltype_{i}" for i in range(ordered.shape[1])]

    # Sanitise cell-type names: h5py treats '/' as a group separator
    ct_names_safe = [n.replace("/", "|") for n in ct_names]

    spatial_adata_full.obsm["cell_type_proportions"] = pd.DataFrame(
        ordered, index=spatial_adata_full.obs_names, columns=ct_names_safe
    )

    # Argmax label (use original names for readability in obs)
    argmax_idx = np.argmax(ordered, axis=1)
    spatial_adata_full.obs[f"{method}_argmax"] = [ct_names[i] for i in argmax_idx]
    log.info(f"Stored obsm['cell_type_proportions'] and obs['{method}_argmax']")

    # --- Save ---
    if output_file is not None:
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        spatial_adata_full.write_h5ad(output_file)
        log.info(f"Saved deconvolution result to {output_file}")

    return spatial_adata_full


# ---------------------------------------------------------------------------
# Existing helper: select_signature_genes (unchanged)
# ---------------------------------------------------------------------------


def select_signature_genes(
    sc_adata: ad.AnnData,
    celltype_key: str,
    sc_batch_key: str,
    n_genes_max: int,
    skip_hvg: bool = True,
    cache_dir: str | None = None,
    force_recompute: bool = False,
    sc_data_file: str | None = None,
    logger_instance: logging.Logger | None = None,
) -> list[str]:
    """Select signature genes for deconvolution.

    This function selects genes for deconvolution by combining:
    1. Highly variable genes (HVGs) - optional
    2. Top marker genes per cell type (always computed)

    Results can be cached to avoid recomputation.

    Parameters
    ----------
    sc_adata : AnnData
        Single-cell reference data
    celltype_key : str
        Column name in sc_adata.obs containing cell type annotations
    sc_batch_key : str
        Column name in sc_adata.obs containing batch information
    n_genes_max : int
        Maximum number of genes to return
    skip_hvg : bool
        If True, skip HVG computation and use only marker genes (faster, less memory)
    cache_dir : str, optional
        Directory to cache signature genes. If None, caching is disabled.
    force_recompute : bool
        If True, force recomputation even if cache exists.
    sc_data_file : str, optional
        Path to single-cell data file (for cache key generation).
    logger_instance : Logger, optional
        Custom logger instance. If None, uses module logger.

    Returns
    -------
    list of str
        List of signature gene names
    """
    log = logger_instance if logger_instance is not None else logger
    log.info("Selecting signature genes...")

    # Check cache if enabled
    if cache_dir and not force_recompute and sc_data_file:
        cache_key = get_cache_key(sc_data_file, celltype_key, sc_batch_key, n_genes_max, skip_hvg)
        cache_path = Path(cache_dir) / f"signature_genes_{cache_key}.pkl"

        cached_genes = load_cached_signatures(cache_path, logger_instance=log)
        if cached_genes is not None:
            return cached_genes
        log.info("   Cache miss or invalid, computing signature genes...")
    elif force_recompute:
        log.info("   Force recompute enabled, computing signature genes...")

    candidate_genes = set()

    # Step 1: Check if HVGs are already computed
    if not skip_hvg and "highly_variable" in sc_adata.var.columns:
        log.info("   Using pre-computed highly variable genes...")
        hvg_genes = set(sc_adata.var[sc_adata.var.highly_variable].index)
        candidate_genes.update(hvg_genes)
        log.info(f"   Found {len(hvg_genes)} pre-computed HVGs")
    elif not skip_hvg:
        # Try to compute HVGs with error handling
        try:
            log.info("   Computing highly variable genes...")
            # Use a simpler method that's less memory intensive
            sc.pp.highly_variable_genes(
                sc_adata, flavor="seurat", n_top_genes=2000, subset=False, batch_key=sc_batch_key
            )
            hvg_genes = set(sc_adata.var[sc_adata.var.highly_variable].index)
            candidate_genes.update(hvg_genes)
            log.info(f"   Found {len(hvg_genes)} HVGs")
        except Exception as e:
            log.warning(f"   Failed to compute HVGs: {e}")
            log.warning("   Continuing with marker genes only...")

    # Step 2: Top marker genes per cell type (always compute these)
    log.info("   Computing marker genes per cell type...")
    try:
        sc.tl.rank_genes_groups(sc_adata, groupby=celltype_key, method="wilcoxon", use_raw=False)
        marker_genes = []
        for group in sc_adata.obs[celltype_key].unique():
            df = sc.get.rank_genes_groups_df(sc_adata, group=group)
            top_genes = df.head(100).names.tolist()
            marker_genes.extend(top_genes)
        marker_genes = set(marker_genes)
        candidate_genes.update(marker_genes)
        log.info(f"   Found {len(marker_genes)} marker genes")
    except Exception as e:
        log.error(f"   Failed to compute marker genes: {e}")
        raise

    # Step 3: If no genes selected, use top expressed genes as fallback
    if len(candidate_genes) == 0:
        log.warning("   No genes selected, using top expressed genes as fallback...")
        # Calculate mean expression per gene
        if hasattr(sc_adata.X, "toarray"):
            mean_expr = np.array(sc_adata.X.mean(axis=0)).flatten()
        else:
            mean_expr = np.array(sc_adata.X.mean(axis=0)).flatten()
        top_indices = np.argsort(mean_expr)[-n_genes_max:][::-1]
        candidate_genes = set(sc_adata.var_names[top_indices])
        log.info(f"   Selected {len(candidate_genes)} top expressed genes")

    # Step 4: Limit to max genes
    candidate_genes = list(candidate_genes)
    if len(candidate_genes) > n_genes_max:
        # Prioritize marker genes if we have both
        if "marker_genes" in locals() and len(marker_genes) > 0:
            # Keep all marker genes, then fill with HVGs or top expressed
            marker_list = list(marker_genes)
            other_genes = [g for g in candidate_genes if g not in marker_genes]
            candidate_genes = marker_list + other_genes[: n_genes_max - len(marker_list)]
        else:
            candidate_genes = candidate_genes[:n_genes_max]

    log.info(f"   Using {len(candidate_genes)} signature genes")

    # Save to cache if enabled
    if cache_dir and sc_data_file:
        cache_key = get_cache_key(sc_data_file, celltype_key, sc_batch_key, n_genes_max, skip_hvg)
        cache_path = Path(cache_dir) / f"signature_genes_{cache_key}.pkl"
        metadata = {
            "sc_data_file": sc_data_file,
            "celltype_key": celltype_key,
            "sc_batch_key": sc_batch_key,
            "n_genes_max": n_genes_max,
            "skip_hvg": skip_hvg,
            "n_genes": len(candidate_genes),
        }
        save_cached_signatures(candidate_genes, cache_path, metadata, logger_instance=log)

    return candidate_genes
