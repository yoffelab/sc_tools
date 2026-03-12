#!/usr/bin/env python3
"""Milestone 1: CosMx 1k vs Xenium MT cross-platform integration benchmark.

Scope: 16 CosMx 1k + 16 Xenium MT samples (same 16 patients, Ileum+Rectum).
This is the KEY milestone — real cross-platform batch effects + CD+UC biology.

Usage:
    python run_m1_benchmark.py

Reads from: data/cosmx_1k_*/adata.p0.h5ad + data/xenium_mt_*/adata.p0.h5ad
Writes to:  results/m1_benchmark/
"""

import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

WORKDIR = Path("/home/fs01/juk4007/elementolab/projects/ibd_spatial")
OUTDIR = WORKDIR / "results" / "m1_benchmark"
FIGDIR = WORKDIR / "figures" / "m1"

OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)

# --- 1. Load and concatenate ---
print("=== Loading samples ===")
adatas = []
for panel, prefix, n_samples in [("cosmx_1k", "cosmx_1k", 16), ("xenium_mt", "xenium_mt", 16)]:
    for i in range(1, n_samples + 1):
        sample_id = f"{prefix}_{i:02d}"
        path = WORKDIR / "data" / sample_id / "adata.p0.h5ad"
        if not path.exists():
            print(f"  MISSING: {path}")
            continue
        a = ad.read_h5ad(path)
        a.obs["platform"] = "CosMx" if "cosmx" in panel else "Xenium"
        a.obs["panel"] = panel
        a.obs["sample"] = sample_id
        a.var_names_make_unique()
        adatas.append(a)
        print(f"  {sample_id}: {a.shape[0]} cells x {a.shape[1]} genes")

print(f"\nLoaded {len(adatas)} samples")

# --- 2. Find shared genes and concatenate ---
print("\n=== Gene intersection ===")
gene_sets = [set(a.var_names) for a in adatas]
shared_genes = gene_sets[0]
for gs in gene_sets[1:]:
    shared_genes = shared_genes & gs
shared_genes = sorted(shared_genes)
print(f"Shared genes across all samples: {len(shared_genes)}")

# Subset to shared genes before concat
adatas_shared = [a[:, shared_genes].copy() for a in adatas]

adata = ad.concat(adatas_shared, join="inner", merge="same")
adata.obs_names_make_unique()
print(f"\nConcatenated: {adata.shape[0]} cells x {adata.shape[1]} genes")
print(f"Platform: {adata.obs['platform'].value_counts().to_dict()}")
print(f"Patients: {adata.obs['patient_id'].nunique()} unique")
if "disease" in adata.obs:
    print(f"Disease: {adata.obs['disease'].value_counts().to_dict()}")
if "tissue_type" in adata.obs:
    print(f"Tissue: {adata.obs['tissue_type'].value_counts().to_dict()}")

# --- 3. QC ---
print("\n=== QC ===")
n_genes = adata.n_vars
pct_top = [t for t in [50, 100, 200] if t < n_genes]
sc.pp.calculate_qc_metrics(adata, percent_top=pct_top, inplace=True)
print(f"Median counts/cell: {adata.obs['total_counts'].median():.0f}")
print(f"Median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}")

# Per-platform QC summary
for platform in adata.obs["platform"].unique():
    mask = adata.obs["platform"] == platform
    med_counts = adata.obs.loc[mask, "total_counts"].median()
    med_genes = adata.obs.loc[mask, "n_genes_by_counts"].median()
    print(f"  {platform}: median_counts={med_counts:.0f}, median_genes={med_genes:.0f}, n_cells={mask.sum()}")

# Filter
min_counts = 10
min_genes = 5
n_before = adata.n_obs
sc.pp.filter_cells(adata, min_counts=min_counts)
sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=10)
print(f"After filtering: {adata.n_obs} cells (removed {n_before - adata.n_obs}), {adata.n_vars} genes")

# --- 4. Preprocessing ---
print("\n=== Preprocessing ===")
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)
n_comps = min(50, adata.n_vars - 1)
sc.tl.pca(adata, n_comps=n_comps)
print(f"PCA: {adata.obsm['X_pca'].shape}")

# Unintegrated baseline
sc.pp.neighbors(adata, use_rep="X_pca", key_added="neighbors_pca")
sc.tl.umap(adata, neighbors_key="neighbors_pca")
adata.obsm["X_umap_pca"] = adata.obsm["X_umap"].copy()
sc.tl.leiden(adata, neighbors_key="neighbors_pca", key_added="leiden_pca", resolution=0.5)

# --- 5. Integration methods ---
print("\n=== Integration ===")
results = {}

# 5a. Unintegrated PCA (baseline)
results["PCA"] = {"embedding_key": "X_pca", "neighbors_key": "neighbors_pca"}

# 5b. Harmony
try:
    import harmonypy

    t0 = time.time()
    adata_harmony = adata.copy()
    ho = harmonypy.run_harmony(adata_harmony.obsm["X_pca"], adata_harmony.obs, "platform")
    emb = ho.Z_corr
    if hasattr(emb, "detach"):
        emb = emb.detach().cpu().numpy()
    if emb.shape[0] != adata_harmony.n_obs:
        emb = emb.T
    adata.obsm["X_harmony"] = np.ascontiguousarray(emb, dtype=np.float32)
    sc.pp.neighbors(adata, use_rep="X_harmony", key_added="neighbors_harmony")
    sc.tl.umap(adata, neighbors_key="neighbors_harmony")
    adata.obsm["X_umap_harmony"] = adata.obsm["X_umap"].copy()
    sc.tl.leiden(adata, neighbors_key="neighbors_harmony", key_added="leiden_harmony", resolution=0.5)
    results["Harmony"] = {"embedding_key": "X_harmony", "neighbors_key": "neighbors_harmony"}
    print(f"  Harmony: done ({time.time()-t0:.0f}s)")
except Exception as e:
    print(f"  Harmony: FAILED ({e})")

# 5c. scVI (requires raw counts)
try:
    import scvi

    t0 = time.time()
    adata_scvi = adata.copy()
    adata_scvi.X = adata_scvi.layers["counts"].copy()
    scvi.model.SCVI.setup_anndata(adata_scvi, batch_key="platform")
    model = scvi.model.SCVI(adata_scvi, n_latent=10, n_hidden=64, n_layers=1)
    # Use early stopping with validation split
    model.train(
        max_epochs=50,
        early_stopping=True,
        early_stopping_patience=5,
        batch_size=512,
        train_size=0.9,
    )
    latent = model.get_latent_representation()
    adata.obsm["X_scvi"] = latent.astype(np.float32)
    sc.pp.neighbors(adata, use_rep="X_scvi", key_added="neighbors_scvi")
    sc.tl.umap(adata, neighbors_key="neighbors_scvi")
    adata.obsm["X_umap_scvi"] = adata.obsm["X_umap"].copy()
    sc.tl.leiden(adata, neighbors_key="neighbors_scvi", key_added="leiden_scvi", resolution=0.5)
    results["scVI"] = {"embedding_key": "X_scvi", "neighbors_key": "neighbors_scvi"}
    print(f"  scVI: done ({time.time()-t0:.0f}s)")
except Exception as e:
    print(f"  scVI: FAILED ({e})")

# 5d. BBKNN
try:
    import bbknn

    t0 = time.time()
    adata_bbknn = adata.copy()
    bbknn.bbknn(adata_bbknn, batch_key="platform", use_rep="X_pca")
    adata.obsp["bbknn_connectivities"] = adata_bbknn.obsp["connectivities"].copy()
    adata.obsp["bbknn_distances"] = adata_bbknn.obsp["distances"].copy()
    adata.uns["neighbors_bbknn"] = {
        "connectivities_key": "bbknn_connectivities",
        "distances_key": "bbknn_distances",
        "params": {"use_rep": "X_pca", "method": "bbknn"},
    }
    sc.tl.umap(adata, neighbors_key="neighbors_bbknn")
    adata.obsm["X_umap_bbknn"] = adata.obsm["X_umap"].copy()
    sc.tl.leiden(adata, neighbors_key="neighbors_bbknn", key_added="leiden_bbknn", resolution=0.5)
    results["BBKNN"] = {"embedding_key": "X_pca", "neighbors_key": "neighbors_bbknn"}
    print(f"  BBKNN: done ({time.time()-t0:.0f}s)")
except Exception as e:
    print(f"  BBKNN: FAILED ({e})")

# 5e. Scanorama
try:
    import scanorama

    t0 = time.time()
    # Split by platform
    platforms = sorted(adata.obs["platform"].unique())
    adatas_by_platform = [adata[adata.obs["platform"] == p].copy() for p in platforms]
    # Scanorama needs log-normalized data in .X (already done above before scale)
    # Re-normalize from counts for Scanorama since .X is now scaled
    for a in adatas_by_platform:
        a.X = a.layers["counts"].copy()
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    corrected, _ = scanorama.correct_scanpy(adatas_by_platform)
    # Reassemble in original order
    corrected_dict = {}
    for p, a in zip(platforms, corrected):
        for idx in a.obs_names:
            corrected_dict[idx] = a[idx].obsm["X_scanorama"][0]
    corrected_X = np.vstack([corrected_dict[idx] for idx in adata.obs_names])
    adata.obsm["X_scanorama"] = corrected_X.astype(np.float32)
    sc.pp.neighbors(adata, use_rep="X_scanorama", key_added="neighbors_scanorama")
    sc.tl.umap(adata, neighbors_key="neighbors_scanorama")
    adata.obsm["X_umap_scanorama"] = adata.obsm["X_umap"].copy()
    sc.tl.leiden(adata, neighbors_key="neighbors_scanorama", key_added="leiden_scanorama", resolution=0.5)
    results["Scanorama"] = {"embedding_key": "X_scanorama", "neighbors_key": "neighbors_scanorama"}
    print(f"  Scanorama: done ({time.time()-t0:.0f}s)")
except Exception as e:
    print(f"  Scanorama: FAILED ({e})")

# --- 6. Benchmark metrics ---
print("\n=== Benchmark ===")
from sklearn.metrics import silhouette_score


def compute_batch_asw(adata, embedding_key, batch_key="platform"):
    """ASW for batch: lower = better mixing."""
    emb = adata.obsm[embedding_key]
    labels = adata.obs[batch_key].values
    if len(np.unique(labels)) < 2:
        return np.nan
    asw = silhouette_score(emb, labels, sample_size=min(10000, adata.n_obs), random_state=42)
    return asw


def compute_patient_mixing(adata, embedding_key, patient_key="patient_id", batch_key="platform"):
    """Per-patient cross-platform mixing score.
    For matched design: within each patient, platforms should mix (ASW ~ 0).
    """
    scores = []
    for patient in adata.obs[patient_key].unique():
        mask = adata.obs[patient_key] == patient
        sub = adata[mask]
        if len(sub.obs[batch_key].unique()) < 2:
            continue
        if sub.n_obs < 20:
            continue
        emb = sub.obsm[embedding_key]
        asw = silhouette_score(emb, sub.obs[batch_key].values, random_state=42)
        scores.append(asw)
    return np.mean(scores) if scores else np.nan


def compute_celltype_asw(adata, embedding_key, celltype_key="celltype"):
    """ASW for cell types: higher = better biological preservation."""
    if celltype_key not in adata.obs:
        return np.nan
    labels = adata.obs[celltype_key].values
    valid = ~pd.isna(labels)
    if valid.sum() < 100:
        return np.nan
    unique_labels = np.unique(labels[valid])
    if len(unique_labels) < 2:
        return np.nan
    emb = adata[valid].obsm[embedding_key]
    asw = silhouette_score(emb, labels[valid], sample_size=min(10000, valid.sum()), random_state=42)
    return asw


benchmark_rows = []
for method, info in results.items():
    emb_key = info["embedding_key"]
    batch_asw = compute_batch_asw(adata, emb_key)
    patient_mix = compute_patient_mixing(adata, emb_key)
    ct_asw = compute_celltype_asw(adata, emb_key, "celltype")
    ct_broad_asw = compute_celltype_asw(adata, emb_key, "celltype_broad")

    row = {
        "method": method,
        "batch_asw": batch_asw,
        "patient_mixing_asw": patient_mix,
        "celltype_asw": ct_asw,
        "celltype_broad_asw": ct_broad_asw,
        "batch_score": 1 - abs(batch_asw),  # higher = better mixing
    }
    benchmark_rows.append(row)
    print(
        f"  {method}: batch_ASW={batch_asw:.3f}, patient_mix={patient_mix:.3f}, "
        f"ct_ASW={ct_asw:.3f}, ct_broad_ASW={ct_broad_asw:.3f}, "
        f"batch_score={row['batch_score']:.3f}"
    )

benchmark_df = pd.DataFrame(benchmark_rows).sort_values("batch_score", ascending=False)
benchmark_df.to_csv(OUTDIR / "m1_benchmark.csv", index=False)
print(f"\nBenchmark saved to {OUTDIR / 'm1_benchmark.csv'}")
print(benchmark_df.to_string(index=False))

# --- 7. Disease composition analysis ---
print("\n=== Disease composition ===")
if "disease" in adata.obs and "celltype" in adata.obs:
    # Cell type proportions by disease
    ct_disease = pd.crosstab(adata.obs["disease"], adata.obs["celltype"], normalize="index")
    ct_disease.to_csv(OUTDIR / "m1_celltype_by_disease.csv")
    print("Cell type proportions by disease:")
    print(ct_disease.round(3).to_string())

    # Fisher test: CD vs UC cell type proportions
    from scipy.stats import fisher_exact

    for ct in adata.obs["celltype"].unique():
        cd_mask = adata.obs["disease"] == "CD"
        uc_mask = adata.obs["disease"] == "UC"
        ct_mask = adata.obs["celltype"] == ct
        table = [
            [(cd_mask & ct_mask).sum(), (cd_mask & ~ct_mask).sum()],
            [(uc_mask & ct_mask).sum(), (uc_mask & ~ct_mask).sum()],
        ]
        if min(table[0][0], table[1][0]) > 0:
            _, p = fisher_exact(table)
            if p < 0.05:
                cd_frac = table[0][0] / (table[0][0] + table[0][1])
                uc_frac = table[1][0] / (table[1][0] + table[1][1])
                print(f"  {ct}: CD={cd_frac:.3f} vs UC={uc_frac:.3f} (Fisher p={p:.4f})")

# --- 8. UMAP plots ---
print("\n=== Plotting ===")
color_keys = ["platform", "patient_id"]
if "disease" in adata.obs:
    color_keys.append("disease")
if "tissue_type" in adata.obs:
    color_keys.append("tissue_type")
if "celltype" in adata.obs:
    color_keys.append("celltype")
if "celltype_broad" in adata.obs:
    color_keys.append("celltype_broad")

n_colors = len(color_keys)

for method, info in results.items():
    umap_key = f"X_umap_{method.lower()}"
    if umap_key not in adata.obsm:
        umap_key = "X_umap_pca"

    fig, axes = plt.subplots(1, n_colors, figsize=(6 * n_colors, 5))
    if n_colors == 1:
        axes = [axes]
    fig.suptitle(f"M1: {method}", fontsize=14, fontweight="bold")

    adata.obsm["X_umap"] = adata.obsm[umap_key]
    for j, key in enumerate(color_keys):
        sc.pl.umap(adata, color=key, ax=axes[j], show=False, title=key)

    plt.tight_layout()
    fig.savefig(FIGDIR / f"m1_umap_{method.lower()}.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: m1_umap_{method.lower()}.png")

# Summary UMAP grid
n_methods = len(results)
fig, axes = plt.subplots(n_methods, n_colors, figsize=(6 * n_colors, 5 * n_methods))
if n_methods == 1:
    axes = axes.reshape(1, -1)

for i, (method, info) in enumerate(results.items()):
    umap_key = f"X_umap_{method.lower()}"
    if umap_key not in adata.obsm:
        umap_key = "X_umap_pca"
    adata.obsm["X_umap"] = adata.obsm[umap_key]

    for j, key in enumerate(color_keys):
        sc.pl.umap(adata, color=key, ax=axes[i, j], show=False, title=f"{method} - {key}")

plt.tight_layout()
fig.savefig(FIGDIR / "m1_umap_grid.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved: m1_umap_grid.png")

# --- 9. Save ---
adata.write_h5ad(OUTDIR / "adata.m1.h5ad")
print(f"\nSaved: {OUTDIR / 'adata.m1.h5ad'} ({adata.n_obs} cells x {adata.n_vars} genes)")

# --- 10. Summary ---
print("\n" + "=" * 60)
print("MILESTONE 1 COMPLETE")
print("=" * 60)
best = benchmark_df.iloc[0]
print(f"Best method: {best['method']} (batch_score={best['batch_score']:.3f})")
print(f"M0 reference: batch_score > 0.99 (technical replicates)")
print(f"M1 target: batch_score > 0.65 (cross-platform)")
if best["batch_score"] < 0.65:
    print("WARNING: batch_score < 0.65 -- cross-platform integration may need tuning")
print(f"\nOutputs:")
print(f"  Benchmark: {OUTDIR / 'm1_benchmark.csv'}")
print(f"  Figures: {FIGDIR}")
print(f"  AnnData: {OUTDIR / 'adata.m1.h5ad'}")
