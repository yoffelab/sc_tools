#!/usr/bin/env python3
"""Milestone 2: CosMx 6k vs Xenium 5K cross-platform integration benchmark.

Scope: 4 CosMx 6k + 4 Xenium 5K samples (same 4 patients, Rectum, UC+Healthy).
High-plex panels — expected ~1500-2000 shared genes (vs 119 in M1).
This should significantly improve celltype separation.

Usage:
    python run_m2_benchmark.py

Reads from: data/cosmx_6k_*/adata.p0.h5ad + data/xenium_5k_*/adata.p0.h5ad
Writes to:  results/m2_benchmark/
"""

import sys
import os
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
import matplotlib.pyplot as plt

WORKDIR = Path("/home/fs01/juk4007/elementolab/projects/ibd_spatial")
OUTDIR = WORKDIR / "results" / "m2_benchmark"
FIGDIR = WORKDIR / "figures" / "m2"

OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)

# --- 1. Load and concatenate ---
print("=== Loading samples ===")
adatas = []
for panel, prefix, n_samples in [("cosmx_6k", "cosmx_6k", 4), ("xenium_5k", "xenium_5k", 4)]:
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
if len(adatas) == 0:
    print("ERROR: No samples found. Exiting.")
    sys.exit(1)

# --- 2. Find shared genes and concatenate ---
print("\n=== Gene intersection ===")
gene_sets = [set(a.var_names) for a in adatas]
cosmx_genes = set()
xenium_genes = set()
for a in adatas:
    if a.obs["platform"].iloc[0] == "CosMx":
        cosmx_genes |= set(a.var_names)
    else:
        xenium_genes |= set(a.var_names)

print(f"CosMx 6k unique genes: {len(cosmx_genes)}")
print(f"Xenium 5K unique genes: {len(xenium_genes)}")

shared_genes = sorted(cosmx_genes & xenium_genes)
print(f"Shared genes: {len(shared_genes)}")

# Also compute per-sample intersection (more conservative)
all_shared = gene_sets[0]
for gs in gene_sets[1:]:
    all_shared = all_shared & gs
all_shared = sorted(all_shared)
print(f"Genes shared across ALL samples: {len(all_shared)}")

# Use all-sample intersection for integration
shared_genes = all_shared

# Subset to shared genes before concat
adatas_shared = [a[:, shared_genes].copy() for a in adatas]

adata = ad.concat(adatas_shared, join="inner", merge="same")
adata.obs_names_make_unique()
print(f"\nConcatenated: {adata.shape[0]} cells x {adata.shape[1]} genes")
print(f"Platform: {adata.obs['platform'].value_counts().to_dict()}")
if "patient_id" in adata.obs:
    print(f"Patients: {adata.obs['patient_id'].nunique()} unique: {sorted(adata.obs['patient_id'].unique())}")
if "disease" in adata.obs:
    print(f"Disease: {adata.obs['disease'].value_counts().to_dict()}")

# --- 3. QC ---
print("\n=== QC ===")
n_genes = adata.n_vars
pct_top = [t for t in [50, 100, 200, 500] if t < n_genes]
sc.pp.calculate_qc_metrics(adata, percent_top=pct_top, inplace=True)
print(f"Median counts/cell: {adata.obs['total_counts'].median():.0f}")
print(f"Median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}")

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

# HVG selection for high-plex data (>500 genes)
n_hvg = min(2000, adata.n_vars)
if adata.n_vars > 500:
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, batch_key="platform")
    n_hvg_selected = adata.var["highly_variable"].sum()
    print(f"HVG selection: {n_hvg_selected} / {adata.n_vars} genes")
    # Keep all genes but mark HVG
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
else:
    adata_hvg = adata.copy()
    print(f"All {adata.n_vars} genes used (< 500, no HVG filtering)")

sc.pp.scale(adata_hvg, max_value=10)
n_comps = min(50, adata_hvg.n_vars - 1)
sc.tl.pca(adata_hvg, n_comps=n_comps)
print(f"PCA: {adata_hvg.obsm['X_pca'].shape}")

# Copy PCA back to full adata
adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]

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
    ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, "platform")
    emb = ho.Z_corr
    if hasattr(emb, "detach"):
        emb = emb.detach().cpu().numpy()
    if emb.shape[0] != adata.n_obs:
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

    # For high-plex: use HVG subset for scVI
    if adata.n_vars > 500 and "highly_variable" in adata.var:
        hvg_mask = adata.var["highly_variable"].values
        adata_scvi = adata_scvi[:, hvg_mask].copy()
        print(f"  scVI using {adata_scvi.n_vars} HVG genes")

    scvi.model.SCVI.setup_anndata(adata_scvi, batch_key="platform")
    # Higher capacity for high-plex data
    model = scvi.model.SCVI(adata_scvi, n_latent=20, n_hidden=128, n_layers=2)
    model.train(
        max_epochs=100,
        early_stopping=True,
        early_stopping_patience=10,
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
    print(f"  scVI: done ({time.time()-t0:.0f}s), final loss={model.history['elbo_train'].iloc[-1]:.1f}")
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
    platforms = sorted(adata.obs["platform"].unique())
    adatas_by_platform = [adata[adata.obs["platform"] == p].copy() for p in platforms]
    for a in adatas_by_platform:
        a.X = a.layers["counts"].copy()
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    corrected, _ = scanorama.correct_scanpy(adatas_by_platform)
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
    emb = adata.obsm[embedding_key]
    labels = adata.obs[batch_key].values
    if len(np.unique(labels)) < 2:
        return np.nan
    return silhouette_score(emb, labels, sample_size=min(10000, adata.n_obs), random_state=42)


def compute_patient_mixing(adata, embedding_key, patient_key="patient_id", batch_key="platform"):
    scores = []
    for patient in adata.obs[patient_key].unique():
        mask = adata.obs[patient_key] == patient
        sub = adata[mask]
        if len(sub.obs[batch_key].unique()) < 2 or sub.n_obs < 20:
            continue
        emb = sub.obsm[embedding_key]
        asw = silhouette_score(emb, sub.obs[batch_key].values, random_state=42)
        scores.append(asw)
    return np.mean(scores) if scores else np.nan


def compute_celltype_asw(adata, embedding_key, celltype_key="celltype"):
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
    return silhouette_score(emb, labels[valid], sample_size=min(10000, valid.sum()), random_state=42)


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
        "batch_score": 1 - abs(batch_asw),
        "n_genes": adata.n_vars,
        "n_cells": adata.n_obs,
    }
    benchmark_rows.append(row)
    print(
        f"  {method}: batch_ASW={batch_asw:.3f}, patient_mix={patient_mix:.3f}, "
        f"ct_ASW={ct_asw:.3f}, ct_broad_ASW={ct_broad_asw:.3f}, "
        f"batch_score={row['batch_score']:.3f}"
    )

benchmark_df = pd.DataFrame(benchmark_rows).sort_values("batch_score", ascending=False)
benchmark_df.to_csv(OUTDIR / "m2_benchmark.csv", index=False)
print(f"\nBenchmark saved to {OUTDIR / 'm2_benchmark.csv'}")
print(benchmark_df.to_string(index=False))

# --- 7. M1 vs M2 comparison ---
print("\n=== M1 vs M2 Comparison ===")
m1_bench = WORKDIR / "results" / "m1_benchmark" / "m1_benchmark.csv"
if m1_bench.exists():
    m1_df = pd.read_csv(m1_bench)
    print("  M1 results:")
    for _, row in m1_df.iterrows():
        print(f"    {row['method']}: batch_score={row['batch_score']:.3f}, ct_broad_ASW={row.get('celltype_broad_asw', 'N/A')}")
    print("  M2 results:")
    for _, row in benchmark_df.iterrows():
        print(f"    {row['method']}: batch_score={row['batch_score']:.3f}, ct_broad_ASW={row.get('celltype_broad_asw', 'N/A')}")
    print(f"\n  Gene count: M1=119, M2={adata.n_vars}")
    # Check if celltype ASW improved
    m1_best_ct = m1_df["celltype_broad_asw"].max() if "celltype_broad_asw" in m1_df else None
    m2_best_ct = benchmark_df["celltype_broad_asw"].max()
    if m1_best_ct is not None and not np.isnan(m2_best_ct):
        delta = m2_best_ct - m1_best_ct
        print(f"  Celltype broad ASW improvement: {delta:+.3f} (M1={m1_best_ct:.3f} -> M2={m2_best_ct:.3f})")
else:
    print("  M1 benchmark not found for comparison")

# --- 8. Disease composition ---
print("\n=== Disease composition ===")
if "disease" in adata.obs and "celltype" in adata.obs:
    ct_disease = pd.crosstab(adata.obs["disease"], adata.obs["celltype"], normalize="index")
    ct_disease.to_csv(OUTDIR / "m2_celltype_by_disease.csv")
    print("Cell type proportions by disease:")
    print(ct_disease.round(3).to_string())
elif "disease" in adata.obs:
    print(f"Disease distribution: {adata.obs['disease'].value_counts().to_dict()}")
    print("No celltype column available for disease composition analysis")

# --- 9. Bio evaluation (inline for M2 since smaller dataset) ---
print("\n=== Bio evaluation ===")

# Per-platform celltype ASW
ct_keys = [k for k in ["celltype", "celltype_broad"] if k in adata.obs.columns]
bio_rows = []
for platform in sorted(adata.obs["platform"].unique()):
    mask = adata.obs["platform"] == platform
    sub = adata[mask]
    for ct_key in ct_keys:
        labels = sub.obs[ct_key].values
        valid = ~pd.isna(labels)
        if valid.sum() < 50:
            continue
        unique_labels = np.unique(labels[valid])
        if len(unique_labels) < 2:
            continue
        for method, info in results.items():
            emb_key = info["embedding_key"]
            emb = sub[valid].obsm[emb_key]
            asw = silhouette_score(emb, labels[valid], sample_size=min(5000, valid.sum()), random_state=42)
            bio_rows.append({
                "platform": platform, "celltype_key": ct_key, "method": method,
                "asw": asw, "n_cells": valid.sum(), "n_types": len(unique_labels),
            })
            print(f"  {platform}/{ct_key}/{method}: ASW={asw:.3f} ({len(unique_labels)} types)")

if bio_rows:
    bio_df = pd.DataFrame(bio_rows)
    bio_df.to_csv(OUTDIR / "m2_bio_platform_asw.csv", index=False)

# Cross-platform kNN label transfer
print("\n=== Cross-platform label transfer ===")
from sklearn.neighbors import KNeighborsClassifier

transfer_rows = []
for ct_key in ct_keys:
    labels = adata.obs[ct_key].values
    valid = ~pd.isna(labels)
    if valid.sum() < 50:
        continue

    cosmx_mask = (adata.obs["platform"] == "CosMx").values & valid
    xenium_mask = (adata.obs["platform"] == "Xenium").values & valid

    if cosmx_mask.sum() < 20 or xenium_mask.sum() < 20:
        continue

    for method, info in results.items():
        emb_key = info["embedding_key"]
        emb = adata.obsm[emb_key]

        knn = KNeighborsClassifier(n_neighbors=15, metric="euclidean")
        knn.fit(emb[cosmx_mask], labels[cosmx_mask])
        acc_c2x = knn.score(emb[xenium_mask], labels[xenium_mask])

        knn2 = KNeighborsClassifier(n_neighbors=15, metric="euclidean")
        knn2.fit(emb[xenium_mask], labels[xenium_mask])
        acc_x2c = knn2.score(emb[cosmx_mask], labels[cosmx_mask])

        transfer_rows.append({
            "method": method, "celltype_key": ct_key,
            "acc_cosmx_to_xenium": acc_c2x, "acc_xenium_to_cosmx": acc_x2c,
            "mean_transfer_acc": (acc_c2x + acc_x2c) / 2,
        })
        print(f"  {method}/{ct_key}: CosMx->Xenium={acc_c2x:.3f}, Xenium->CosMx={acc_x2c:.3f}")

if transfer_rows:
    transfer_df = pd.DataFrame(transfer_rows)
    transfer_df.to_csv(OUTDIR / "m2_bio_label_transfer.csv", index=False)

# --- 10. UMAP plots ---
print("\n=== Plotting ===")
color_keys = ["platform", "patient_id"]
if "disease" in adata.obs:
    color_keys.append("disease")
if "celltype" in adata.obs:
    color_keys.append("celltype")
if "celltype_broad" in adata.obs:
    color_keys.append("celltype_broad")

n_colors = len(color_keys)

for method in results:
    umap_key = f"X_umap_{method.lower()}"
    if umap_key not in adata.obsm:
        umap_key = "X_umap_pca"

    fig, axes = plt.subplots(1, n_colors, figsize=(6 * n_colors, 5))
    if n_colors == 1:
        axes = [axes]
    fig.suptitle(f"M2: {method}", fontsize=14, fontweight="bold")

    adata.obsm["X_umap"] = adata.obsm[umap_key]
    for j, key in enumerate(color_keys):
        sc.pl.umap(adata, color=key, ax=axes[j], show=False, title=key)

    plt.tight_layout()
    fig.savefig(FIGDIR / f"m2_umap_{method.lower()}.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: m2_umap_{method.lower()}.png")

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
fig.savefig(FIGDIR / "m2_umap_grid.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved: m2_umap_grid.png")

# --- 11. Save ---
adata.write_h5ad(OUTDIR / "adata.m2.h5ad")
print(f"\nSaved: {OUTDIR / 'adata.m2.h5ad'} ({adata.n_obs} cells x {adata.n_vars} genes)")

# --- 12. Summary ---
print("\n" + "=" * 60)
print("MILESTONE 2 COMPLETE")
print("=" * 60)
best = benchmark_df.iloc[0]
print(f"Best method: {best['method']} (batch_score={best['batch_score']:.3f})")
print(f"Shared genes: {adata.n_vars} (vs M1: 119, M0: 377)")
print(f"Cells: {adata.n_obs}")
if not np.isnan(best["celltype_broad_asw"]):
    print(f"Celltype broad ASW: {best['celltype_broad_asw']:.3f}")
print(f"\nOutputs:")
print(f"  Benchmark: {OUTDIR / 'm2_benchmark.csv'}")
print(f"  Figures: {FIGDIR}")
print(f"  AnnData: {OUTDIR / 'adata.m2.h5ad'}")
