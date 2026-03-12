#!/usr/bin/env python3
"""M2 gimVI: gene imputation across panels using scRNA-seq reference.

gimVI (scvi.model.GIMVI) jointly models spatial + scRNA-seq data to impute
genes not present in one panel. This bridges the CosMx 6k / Xenium 5K gap
by leveraging the SAHA_IBD_RNA.h5ad scRNA-seq reference.

Reads:
  results/m2_benchmark/adata.m2.h5ad
  data/reference/SAHA_IBD_RNA.h5ad
Writes:
  results/m2_benchmark/adata.m2.gimvi.h5ad
  results/m2_benchmark/m2_gimvi_metrics.csv
"""

import os
import sys
import time
from pathlib import Path

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout.reconfigure(line_buffering=True)

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib

matplotlib.use("Agg")

WORKDIR = Path("/home/fs01/juk4007/elementolab/projects/ibd_spatial")
OUTDIR = WORKDIR / "results" / "m2_benchmark"
FIGDIR = WORKDIR / "figures" / "m2"

# --- 1. Load spatial data ---
print("=== Loading M2 spatial data ===")
adata_spatial = ad.read_h5ad(OUTDIR / "adata.m2.h5ad")
print(f"Spatial: {adata_spatial.n_obs} cells x {adata_spatial.n_vars} genes")

# --- 2. Load scRNA-seq reference ---
print("\n=== Loading scRNA-seq reference ===")
ref_path = WORKDIR / "data" / "reference" / "SAHA_IBD_RNA.h5ad"
if not ref_path.exists():
    print(f"ERROR: Reference not found at {ref_path}")
    # Try alternate path
    ref_path = Path("/athena/project-saha/data_IBD/SAHA_IBD_RNA.h5ad")
    if not ref_path.exists():
        print(f"ERROR: Reference not found at {ref_path} either. Exiting.")
        sys.exit(1)

adata_ref = ad.read_h5ad(ref_path)
print(f"Reference: {adata_ref.n_obs} cells x {adata_ref.n_vars} genes")
print(f"Reference obs columns: {list(adata_ref.obs.columns[:20])}")

# Check for celltype columns
for col in ["ct_major_new", "ct_minor_new", "celltype", "cell_type"]:
    if col in adata_ref.obs:
        print(f"  {col}: {adata_ref.obs[col].nunique()} types")

# --- 3. Find shared genes ---
print("\n=== Gene intersection ===")
shared_genes = sorted(set(adata_spatial.var_names) & set(adata_ref.var_names))
print(f"Shared genes (spatial x reference): {len(shared_genes)}")
print(f"  Spatial-only: {adata_spatial.n_vars - len(shared_genes)}")
print(f"  Reference-only: {adata_ref.n_vars - len(shared_genes)}")

if len(shared_genes) < 100:
    print("ERROR: Too few shared genes for gimVI. Exiting.")
    sys.exit(1)

# Subset both to shared genes
adata_spatial_sub = adata_spatial[:, shared_genes].copy()
adata_ref_sub = adata_ref[:, shared_genes].copy()

# Ensure raw counts
if "counts" in adata_spatial_sub.layers:
    adata_spatial_sub.X = adata_spatial_sub.layers["counts"].copy()
    print("  Using counts layer for spatial data")

# Check if reference has raw counts
x_max = adata_ref_sub.X.max() if hasattr(adata_ref_sub.X, 'max') else np.max(adata_ref_sub.X.toarray()[:100])
print(f"  Reference X max (sample): {x_max:.1f}")
if x_max < 20:
    print("  WARNING: Reference may be normalized, not raw counts. gimVI prefers raw counts.")
    if "counts" in adata_ref_sub.layers:
        adata_ref_sub.X = adata_ref_sub.layers["counts"].copy()
        print("  Using counts layer for reference")
    elif adata_ref_sub.raw is not None:
        print("  Trying adata.raw for reference counts")
        adata_ref_raw = adata_ref_sub.raw.to_adata()
        shared_in_raw = sorted(set(shared_genes) & set(adata_ref_raw.var_names))
        if len(shared_in_raw) > len(shared_genes) * 0.8:
            adata_ref_sub = adata_ref_raw[:, shared_in_raw].copy()
            adata_spatial_sub = adata_spatial_sub[:, shared_in_raw].copy()
            print(f"  Using raw reference: {len(shared_in_raw)} shared genes")

# Subsample reference if too large
MAX_REF_CELLS = 50000
if adata_ref_sub.n_obs > MAX_REF_CELLS:
    np.random.seed(42)
    idx = np.random.choice(adata_ref_sub.n_obs, MAX_REF_CELLS, replace=False)
    adata_ref_sub = adata_ref_sub[idx].copy()
    print(f"  Subsampled reference to {MAX_REF_CELLS} cells")

print(f"\nFinal: spatial={adata_spatial_sub.shape}, reference={adata_ref_sub.shape}")

# --- 4. Run gimVI ---
print("\n=== Running gimVI ===")
try:
    import scvi
    from scvi.model import GIMVI

    t0 = time.time()

    GIMVI.setup_anndata(adata_spatial_sub)
    GIMVI.setup_anndata(adata_ref_sub)

    model = GIMVI(adata_spatial_sub, adata_ref_sub, n_latent=20)
    model.train(max_epochs=100, batch_size=256)

    # Get latent representations
    latent_spatial, latent_ref = model.get_latent_representation()
    print(f"  gimVI latent spatial: {latent_spatial.shape}")
    print(f"  gimVI latent reference: {latent_ref.shape}")

    # Get imputed values for spatial data
    imputed = model.get_imputed_values(normalized=True)
    print(f"  Imputed expression matrix: {imputed.shape}")

    # Store in spatial adata
    adata_spatial.obsm["X_gimvi"] = latent_spatial.astype(np.float32)

    # Build neighbors and UMAP from gimVI latent
    sc.pp.neighbors(adata_spatial, use_rep="X_gimvi", key_added="neighbors_gimvi")
    sc.tl.umap(adata_spatial, neighbors_key="neighbors_gimvi")
    adata_spatial.obsm["X_umap_gimvi"] = adata_spatial.obsm["X_umap"].copy()
    sc.tl.leiden(adata_spatial, neighbors_key="neighbors_gimvi", key_added="leiden_gimvi", resolution=0.5)

    print(f"  gimVI: done ({time.time()-t0:.0f}s)")

    # --- 5. Evaluate ---
    print("\n=== gimVI evaluation ===")
    from sklearn.metrics import silhouette_score

    emb = adata_spatial.obsm["X_gimvi"]

    # Batch ASW
    batch_asw = silhouette_score(emb, adata_spatial.obs["platform"].values,
                                  sample_size=min(10000, adata_spatial.n_obs), random_state=42)
    print(f"  Batch ASW: {batch_asw:.3f} (batch_score={1-abs(batch_asw):.3f})")

    # Celltype ASW
    for ct_key in ["celltype", "celltype_broad"]:
        if ct_key in adata_spatial.obs:
            labels = adata_spatial.obs[ct_key].values
            valid = ~pd.isna(labels)
            if valid.sum() > 100 and len(np.unique(labels[valid])) > 1:
                ct_asw = silhouette_score(emb[valid], labels[valid],
                                          sample_size=min(10000, valid.sum()), random_state=42)
                print(f"  {ct_key} ASW: {ct_asw:.3f}")

    # Platform entropy
    from scipy.stats import entropy as scipy_entropy
    n_platforms = adata_spatial.obs["platform"].nunique()
    max_ent = np.log(n_platforms)
    entropies = []
    for cluster in adata_spatial.obs["leiden_gimvi"].unique():
        mask = adata_spatial.obs["leiden_gimvi"] == cluster
        counts = adata_spatial.obs.loc[mask, "platform"].value_counts()
        ent = scipy_entropy(counts / counts.sum()) / max_ent if max_ent > 0 else 0
        entropies.append(ent)
    print(f"  Median platform entropy: {np.median(entropies):.3f}")

    # Save
    adata_spatial.write_h5ad(OUTDIR / "adata.m2.gimvi.h5ad")
    print(f"\nSaved: {OUTDIR / 'adata.m2.gimvi.h5ad'}")

    # UMAP plot
    import matplotlib.pyplot as plt
    color_keys = ["platform", "patient_id", "disease", "celltype_broad"]
    color_keys = [k for k in color_keys if k in adata_spatial.obs]
    n_colors = len(color_keys)

    fig, axes = plt.subplots(1, n_colors, figsize=(6 * n_colors, 5))
    if n_colors == 1:
        axes = [axes]
    fig.suptitle("M2: gimVI", fontsize=14, fontweight="bold")
    adata_spatial.obsm["X_umap"] = adata_spatial.obsm["X_umap_gimvi"]
    for j, key in enumerate(color_keys):
        sc.pl.umap(adata_spatial, color=key, ax=axes[j], show=False, title=key)
    plt.tight_layout()
    fig.savefig(FIGDIR / "m2_umap_gimvi.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: m2_umap_gimvi.png")

except Exception as e:
    print(f"  gimVI: FAILED ({e})")
    import traceback
    traceback.print_exc()
