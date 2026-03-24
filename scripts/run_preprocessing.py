#!/usr/bin/env python
"""
run_preprocessing.py -- Phase 2 preprocessing for Visium HD filtered datasets.

Methods
-------
harmony   normalize -> log1p -> scale -> HVGs -> PCA -> Harmony -> UMAP -> Leiden
scvi      raw counts -> HVGs -> scVI (GPU) -> latent -> UMAP -> Leiden
resolvi   raw counts -> HVGs -> resolVI (spatial GPU) -> latent -> UMAP -> Leiden

Usage
-----
python run_preprocessing.py \
    --adata  results/adata.robin_008um.filtered.h5ad \
    --method harmony \
    --output results/preprocessing/robin_008um_harmony.h5ad \
    --batch-key sample \
    --n-hvgs 3000
"""

import argparse
import logging
import sys

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Arrow string coercion (needed when reading h5ad written on newer pandas)
# ---------------------------------------------------------------------------
pd.set_option("future.infer_string", False)
try:
    pd.options.mode.string_storage = "python"
except Exception:
    pass


def _coerce_arrow_strings(adata):
    """Convert Arrow-backed StringDtype columns to plain object dtype."""
    for col in adata.obs.columns:
        if hasattr(adata.obs[col].dtype, "name") and "string" in str(adata.obs[col].dtype).lower():
            adata.obs[col] = adata.obs[col].astype(object)
    for col in adata.var.columns:
        if hasattr(adata.var[col].dtype, "name") and "string" in str(adata.var[col].dtype).lower():
            adata.var[col] = adata.var[col].astype(object)
    if hasattr(adata.obs_names.dtype, "name") and "string" in str(adata.obs_names.dtype).lower():
        adata.obs_names = adata.obs_names.astype(object)
    if hasattr(adata.var_names.dtype, "name") and "string" in str(adata.var_names.dtype).lower():
        adata.var_names = adata.var_names.astype(object)
    return adata


# ---------------------------------------------------------------------------
# Try to import RAPIDS
# ---------------------------------------------------------------------------
try:
    import rapids_singlecell as rsc
    import cupy as cp
    HAS_RAPIDS = True
    log.info("rapids_singlecell available — using GPU backend")
except ImportError:
    HAS_RAPIDS = False
    log.warning("rapids_singlecell not available — falling back to CPU (scanpy)")


def _try_gpu(adata):
    """Move adata to GPU if RAPIDS is available."""
    if HAS_RAPIDS:
        rsc.get.anndata_to_GPU(adata)
    return adata


def _to_cpu(adata):
    """Move adata back to CPU."""
    if HAS_RAPIDS:
        rsc.get.anndata_to_CPU(adata)
    return adata


# ---------------------------------------------------------------------------
# HVG selection (subsample for very large datasets to keep it fast)
# ---------------------------------------------------------------------------
HVG_SUBSAMPLE = 100_000  # subsample N obs for HVG fitting if larger


def select_hvgs(adata, n_hvgs: int, batch_key: str | None = None):
    """Select highly variable genes using normalized log1p data.

    Uses 'seurat' flavor (normalized + log1p) — robust across all scanpy versions.
    Works on a subsample for very large datasets; raw X is never modified.
    """
    import scanpy as sc

    n_obs = adata.n_obs

    def _hvg_on(adata_norm):
        sc.pp.highly_variable_genes(
            adata_norm,
            n_top_genes=n_hvgs,
            batch_key=batch_key,
            flavor="seurat",
        )
        return adata_norm.var["highly_variable"]

    if n_obs > HVG_SUBSAMPLE:
        log.info(f"Subsampling {HVG_SUBSAMPLE:,} / {n_obs:,} obs for HVG selection")
        idx = np.random.default_rng(42).choice(n_obs, HVG_SUBSAMPLE, replace=False)
        adata_sub = adata[idx].copy()
        sc.pp.normalize_total(adata_sub, target_sum=1e4)
        sc.pp.log1p(adata_sub)
        adata.var["highly_variable"] = _hvg_on(adata_sub)
    else:
        adata_norm = adata.copy()
        sc.pp.normalize_total(adata_norm, target_sum=1e4)
        sc.pp.log1p(adata_norm)
        adata.var["highly_variable"] = _hvg_on(adata_norm)


def _ensure_obs_unique(adata):
    if not adata.obs_names.is_unique:
        log.warning("obs_names not unique — making unique")
        adata.obs_names_make_unique()
    return adata


# ---------------------------------------------------------------------------
# Harmony pipeline (GPU if available)
# ---------------------------------------------------------------------------
def run_harmony(adata, batch_key: str, n_hvgs: int, n_pcs: int = 50, resolution: float = 0.5):
    import scanpy as sc

    log.info("=== Harmony pipeline ===")
    _ensure_obs_unique(adata)

    # Save raw counts
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
        log.info("Saved raw counts to layers['counts']")

    # HVG selection
    select_hvgs(adata, n_hvgs=n_hvgs, batch_key=batch_key)
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    log.info(f"Selected {adata_hvg.n_vars:,} HVGs")

    # Normalize -> log1p -> scale (zero_center=False to avoid densifying large sparse matrices)
    sc.pp.normalize_total(adata_hvg, target_sum=1e4)
    sc.pp.log1p(adata_hvg)
    sc.pp.scale(adata_hvg, max_value=10, zero_center=False)

    # PCA + Harmony
    if HAS_RAPIDS:
        _try_gpu(adata_hvg)
        rsc.pp.pca(adata_hvg, n_comps=n_pcs)
        rsc.pp.harmony_integrate(adata_hvg, key=batch_key, max_iter_harmony=30)
        use_rep = "X_pca_harmony"
        rsc.pp.neighbors(adata_hvg, use_rep=use_rep, n_neighbors=15)
        rsc.tl.umap(adata_hvg, min_dist=0.3)
        rsc.tl.leiden(adata_hvg, resolution=resolution, key_added="leiden")
        _to_cpu(adata_hvg)
    else:
        sc.pp.pca(adata_hvg, n_comps=n_pcs)
        sc.external.pp.harmony_integrate(adata_hvg, key=batch_key, max_iter_harmony=30)
        use_rep = "X_pca_harmony"
        sc.pp.neighbors(adata_hvg, use_rep=use_rep, n_neighbors=15)
        sc.tl.umap(adata_hvg, min_dist=0.3)
        sc.tl.leiden(adata_hvg, resolution=resolution, key_added="leiden")

    # Transfer embeddings/clusters back to full adata
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata.obsm["X_pca_harmony"] = adata_hvg.obsm["X_pca_harmony"]
    adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
    adata.obs["leiden"] = adata_hvg.obs["leiden"]
    adata.var["highly_variable"] = adata.var.get("highly_variable", False)

    log.info(f"Harmony done: {adata.n_obs:,} obs, {adata.obs['leiden'].nunique()} clusters")
    return adata


# ---------------------------------------------------------------------------
# scVI pipeline
# ---------------------------------------------------------------------------
def run_scvi(adata, batch_key: str, n_hvgs: int, n_latent: int = 30, max_epochs: int | None = None):
    import scanpy as sc
    import scvi

    log.info("=== scVI pipeline ===")
    _ensure_obs_unique(adata)

    # Save raw counts
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
        log.info("Saved raw counts to layers['counts']")

    # HVG selection on raw counts (seurat_v3 flavor expects raw)
    if "highly_variable" not in adata.var.columns:
        select_hvgs(adata, n_hvgs=n_hvgs, batch_key=batch_key)

    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    # Restore raw counts for scVI
    adata_hvg.X = adata_hvg.layers["counts"].copy()
    log.info(f"Using {adata_hvg.n_vars:,} HVGs for scVI")

    scvi.settings.seed = 42
    if max_epochs is None:
        # Heuristic: fewer epochs for very large datasets
        max_epochs = max(100, min(400, int(20_000 / (adata_hvg.n_obs / 10_000))))
    log.info(f"Training scVI for up to {max_epochs} epochs")

    scvi.model.SCVI.setup_anndata(adata_hvg, layer="counts", batch_key=batch_key)
    model = scvi.model.SCVI(adata_hvg, n_latent=n_latent, n_layers=2, n_hidden=128)
    model.train(max_epochs=max_epochs, accelerator="gpu", devices=1)

    adata_hvg.obsm["X_scVI"] = model.get_latent_representation()
    log.info("scVI training complete; computing UMAP")

    if HAS_RAPIDS:
        _try_gpu(adata_hvg)
        rsc.pp.neighbors(adata_hvg, use_rep="X_scVI", n_neighbors=15)
        rsc.tl.umap(adata_hvg, min_dist=0.3)
        rsc.tl.leiden(adata_hvg, resolution=0.5, key_added="leiden")
        _to_cpu(adata_hvg)
    else:
        sc.pp.neighbors(adata_hvg, use_rep="X_scVI", n_neighbors=15)
        sc.tl.umap(adata_hvg, min_dist=0.3)
        sc.tl.leiden(adata_hvg, resolution=0.5, key_added="leiden")

    # Transfer back to full adata
    adata.obsm["X_scVI"] = adata_hvg.obsm["X_scVI"]
    adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
    adata.obs["leiden"] = adata_hvg.obs["leiden"]
    adata.var["highly_variable"] = adata.var.get("highly_variable", False)

    log.info(f"scVI done: {adata.n_obs:,} obs, {adata.obs['leiden'].nunique()} clusters")
    return adata, model


# ---------------------------------------------------------------------------
# resolVI pipeline
# ---------------------------------------------------------------------------
def run_resolvi(
    adata,
    batch_key: str,
    n_hvgs: int,
    n_latent: int = 20,
    max_epochs: int | None = None,
    spatial_key: str = "spatial",
):
    import scanpy as sc
    import scvi

    log.info("=== resolVI pipeline ===")
    _ensure_obs_unique(adata)

    if spatial_key not in adata.obsm:
        raise ValueError(
            f"resolVI requires spatial coordinates in adata.obsm['{spatial_key}']. "
            "Key not found."
        )

    # Save raw counts
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
        log.info("Saved raw counts to layers['counts']")

    # HVG selection
    if "highly_variable" not in adata.var.columns:
        select_hvgs(adata, n_hvgs=n_hvgs, batch_key=batch_key)

    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    adata_hvg.X = adata_hvg.layers["counts"].copy()

    # resolVI requires all cells to have >0 total counts across HVGs
    import scipy.sparse as sp
    if sp.issparse(adata_hvg.X):
        cell_counts = np.array(adata_hvg.X.sum(axis=1)).ravel()
    else:
        cell_counts = np.array(adata_hvg.X.sum(axis=1)).ravel()
    min_counts = 5  # resolVI recommends >= 5
    keep_mask = cell_counts >= min_counts
    n_removed = (~keep_mask).sum()
    if n_removed > 0:
        log.warning(
            f"Removing {n_removed:,} cells with <{min_counts} total counts across HVGs "
            f"({n_removed / adata_hvg.n_obs * 100:.2f}% of cells) — required by resolVI"
        )
        adata_hvg = adata_hvg[keep_mask].copy()

    log.info(f"Using {adata_hvg.n_vars:,} HVGs for resolVI ({adata_hvg.n_obs:,} cells after filtering)")

    scvi.settings.seed = 42
    if max_epochs is None:
        max_epochs = max(100, min(400, int(20_000 / (adata_hvg.n_obs / 10_000))))
    log.info(f"Training resolVI for up to {max_epochs} epochs")

    # resolVI is in scvi.external for older versions, scvi.model for newer
    try:
        from scvi.external import RESOLVI
        log.info("Using scvi.external.RESOLVI")
    except ImportError:
        from scvi.model import RESOLVI
        log.info("Using scvi.model.RESOLVI")

    # resolVI._prepare_data expects obsm['X_spatial'] (not 'spatial')
    if spatial_key in adata_hvg.obsm and "X_spatial" not in adata_hvg.obsm:
        adata_hvg.obsm["X_spatial"] = adata_hvg.obsm[spatial_key].copy()
        log.info(f"Copied obsm['{spatial_key}'] → obsm['X_spatial'] for resolVI")

    RESOLVI.setup_anndata(
        adata_hvg,
        layer="counts",
        batch_key=batch_key,
    )
    model = RESOLVI(adata_hvg, n_latent=n_latent)
    model.train(max_epochs=max_epochs)

    adata_hvg.obsm["X_resolVI"] = model.get_latent_representation()
    log.info("resolVI training complete; computing UMAP")

    if HAS_RAPIDS:
        _try_gpu(adata_hvg)
        rsc.pp.neighbors(adata_hvg, use_rep="X_resolVI", n_neighbors=15)
        rsc.tl.umap(adata_hvg, min_dist=0.3)
        rsc.tl.leiden(adata_hvg, resolution=0.5, key_added="leiden")
        _to_cpu(adata_hvg)
    else:
        sc.pp.neighbors(adata_hvg, use_rep="X_resolVI", n_neighbors=15)
        sc.tl.umap(adata_hvg, min_dist=0.3)
        sc.tl.leiden(adata_hvg, resolution=0.5, key_added="leiden")

    # Transfer back — adata_hvg may have fewer cells if zero-count cells were removed
    if adata_hvg.n_obs < adata.n_obs:
        # Use index alignment: only cells present in adata_hvg get embeddings
        shared_idx = adata.obs_names.isin(adata_hvg.obs_names)
        n_latent_dim = adata_hvg.obsm["X_resolVI"].shape[1]
        n_umap_dim = adata_hvg.obsm["X_umap"].shape[1]

        resolvi_full = np.full((adata.n_obs, n_latent_dim), np.nan)
        umap_full = np.full((adata.n_obs, n_umap_dim), np.nan)
        resolvi_full[shared_idx] = adata_hvg.obsm["X_resolVI"]
        umap_full[shared_idx] = adata_hvg.obsm["X_umap"]
        adata.obsm["X_resolVI"] = resolvi_full
        adata.obsm["X_umap"] = umap_full

        leiden_full = pd.Series("unassigned", index=adata.obs_names, dtype="object")
        leiden_full.loc[adata_hvg.obs_names] = adata_hvg.obs["leiden"].astype(str).values
        adata.obs["leiden"] = pd.Categorical(leiden_full)

        log.info(
            f"Transferred embeddings for {adata_hvg.n_obs:,} / {adata.n_obs:,} cells "
            f"({(~shared_idx).sum():,} cells marked unassigned due to low counts)"
        )
    else:
        adata.obsm["X_resolVI"] = adata_hvg.obsm["X_resolVI"]
        adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
        adata.obs["leiden"] = adata_hvg.obs["leiden"]
    adata.var["highly_variable"] = adata.var.get("highly_variable", False)

    log.info(f"resolVI done: {adata.n_obs:,} obs, {adata.obs['leiden'].nunique()} clusters")
    return adata, model


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--adata", required=True, help="Path to filtered h5ad input")
    p.add_argument(
        "--method",
        required=True,
        choices=["harmony", "scvi", "resolvi"],
        help="Integration method",
    )
    p.add_argument("--output", required=True, help="Path to write output h5ad")
    p.add_argument("--batch-key", default="sample", help="obs column for batch (default: sample)")
    p.add_argument("--n-hvgs", type=int, default=3000, help="Number of HVGs (default: 3000)")
    p.add_argument("--n-pcs", type=int, default=50, help="PCs for Harmony (default: 50)")
    p.add_argument("--n-latent", type=int, default=30, help="scVI/resolVI latent dims (default: 30)")
    p.add_argument("--max-epochs", type=int, default=None, help="Max training epochs (auto if not set)")
    p.add_argument("--resolution", type=float, default=0.5, help="Leiden resolution (default: 0.5)")
    p.add_argument(
        "--save-model-dir",
        default=None,
        help="Directory to save scVI/resolVI model (optional)",
    )
    return p.parse_args()


def main():
    import pathlib
    import scanpy as sc

    args = parse_args()

    log.info(f"Loading {args.adata}")
    adata = sc.read_h5ad(args.adata)
    adata = _coerce_arrow_strings(adata)
    log.info(f"Loaded: {adata.n_obs:,} obs × {adata.n_vars:,} vars")

    # Validate batch key
    if args.batch_key not in adata.obs.columns:
        available = list(adata.obs.columns)
        raise ValueError(f"batch-key '{args.batch_key}' not in obs. Available: {available}")

    n_batches = adata.obs[args.batch_key].nunique()
    log.info(f"Batch key='{args.batch_key}', {n_batches} unique values")

    if args.method == "harmony":
        adata = run_harmony(
            adata,
            batch_key=args.batch_key,
            n_hvgs=args.n_hvgs,
            n_pcs=args.n_pcs,
            resolution=args.resolution,
        )
        model = None

    elif args.method == "scvi":
        adata, model = run_scvi(
            adata,
            batch_key=args.batch_key,
            n_hvgs=args.n_hvgs,
            n_latent=args.n_latent,
            max_epochs=args.max_epochs,
        )

    elif args.method == "resolvi":
        adata, model = run_resolvi(
            adata,
            batch_key=args.batch_key,
            n_hvgs=args.n_hvgs,
            n_latent=args.n_latent,
            max_epochs=args.max_epochs,
        )

    # Save output
    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    log.info(f"Writing output to {out_path}")
    adata.write_h5ad(out_path)
    log.info("Done.")

    # Optionally save model
    if model is not None and args.save_model_dir:
        model_dir = pathlib.Path(args.save_model_dir)
        model_dir.mkdir(parents=True, exist_ok=True)
        model.save(str(model_dir), overwrite=True)
        log.info(f"Model saved to {model_dir}")


if __name__ == "__main__":
    main()
