#!/usr/bin/env python3
"""M2 follow-up: scANVI, Scanorama fix, IBD marker check, platform entropy.

Addresses gaps from the original plan:
1. scANVI (semi-supervised) — uses celltype_broad labels during training
2. Scanorama — debug and fix the X_scanorama key error
3. IBD marker gene check — are canonical markers in the 2,552-gene intersection?
4. Per-cluster platform entropy — planned success criterion never computed
5. Formal success criteria evaluation from the plan

Reads: results/m2_benchmark/adata.m2.h5ad
Writes: results/m2_benchmark/m2_followup.csv, updated adata
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
import matplotlib.pyplot as plt

WORKDIR = Path("/home/fs01/juk4007/elementolab/projects/ibd_spatial")
OUTDIR = WORKDIR / "results" / "m2_benchmark"
FIGDIR = WORKDIR / "figures" / "m2"

print("=== Loading M2 adata ===")
adata = ad.read_h5ad(OUTDIR / "adata.m2.h5ad")
print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")
print(f"Existing embeddings: {[k for k in adata.obsm if k.startswith('X_')]}")

# ============================================================
# 1. IBD Marker Gene Check
# ============================================================
print("\n=== 1. IBD Marker Gene Check ===")
canonical_markers = ["EPCAM", "CD3E", "CD68", "MS4A1", "MKI67", "FOXP3", "IFNG", "TNF"]
present = [g for g in canonical_markers if g in adata.var_names]
missing = [g for g in canonical_markers if g not in adata.var_names]
print(f"Canonical IBD markers in 2,552-gene intersection: {len(present)}/8")
print(f"  Present: {present}")
print(f"  Missing: {missing}")
print(f"  Criterion (>=5/8): {'PASS' if len(present) >= 5 else 'FAIL'}")

# Also check additional relevant markers
extra_markers = {
    "Epithelial": ["KRT20", "MUC2", "LGR5", "OLFM4", "CDH1"],
    "T cell": ["CD4", "CD8A", "CD8B", "GZMB", "PRF1", "IL17A", "IL22", "RORC"],
    "Myeloid": ["CD14", "CD163", "ITGAX", "CSF1R", "IL1B", "IL6"],
    "B/Plasma": ["CD19", "CD79A", "JCHAIN", "MZB1", "SDC1"],
    "Fibroblast": ["COL1A1", "ACTA2", "VIM", "FAP", "THY1"],
    "IBD-specific": ["NOD2", "ATG16L1", "IL23R", "CARD9", "REG1A", "DMBT1", "LCN2"],
}
print("\nAdditional marker presence:")
for category, markers in extra_markers.items():
    found = [g for g in markers if g in adata.var_names]
    print(f"  {category}: {len(found)}/{len(markers)} — {found}")

# ============================================================
# 2. scANVI (semi-supervised scVI)
# ============================================================
print("\n=== 2. scANVI ===")
try:
    import scvi

    t0 = time.time()

    # Prepare data — scANVI needs raw counts and celltype labels
    adata_scanvi = adata.copy()
    adata_scanvi.X = adata_scanvi.layers["counts"].copy()

    # Use HVG subset like scVI did
    if "highly_variable" in adata.var:
        hvg_mask = adata.var["highly_variable"].values
        adata_scanvi = adata_scanvi[:, hvg_mask].copy()
        print(f"  Using {adata_scanvi.n_vars} HVG genes")

    # Check label coverage
    ct_key = "celltype_broad"
    labels = adata_scanvi.obs[ct_key].values
    valid = ~pd.isna(labels)
    print(f"  Label coverage: {valid.sum()}/{adata_scanvi.n_obs} cells ({100*valid.sum()/adata_scanvi.n_obs:.1f}%)")

    # Handle NaN labels — scANVI needs a special "unlabeled" category
    # Convert from Categorical to str first to avoid "new category" error
    adata_scanvi.obs["_scanvi_labels"] = adata_scanvi.obs[ct_key].astype(str).replace("nan", "Unknown")

    scvi.model.SCVI.setup_anndata(adata_scanvi, batch_key="platform")
    scvi_model = scvi.model.SCVI(adata_scanvi, n_latent=20, n_hidden=128, n_layers=2)
    scvi_model.train(
        max_epochs=50,
        early_stopping=True,
        early_stopping_patience=10,
        batch_size=512,
        train_size=0.9,
    )
    print(f"  scVI pretrain done ({time.time()-t0:.0f}s)")

    # Initialize scANVI from pretrained scVI
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        unlabeled_category="Unknown",
        labels_key="_scanvi_labels",
    )
    scanvi_model.train(
        max_epochs=30,
        early_stopping=True,
        early_stopping_patience=5,
        batch_size=512,
        train_size=0.9,
    )
    latent = scanvi_model.get_latent_representation()
    adata.obsm["X_scanvi"] = latent.astype(np.float32)

    # Neighbors, UMAP, leiden
    sc.pp.neighbors(adata, use_rep="X_scanvi", key_added="neighbors_scanvi")
    sc.tl.umap(adata, neighbors_key="neighbors_scanvi")
    adata.obsm["X_umap_scanvi"] = adata.obsm["X_umap"].copy()
    sc.tl.leiden(adata, neighbors_key="neighbors_scanvi", key_added="leiden_scanvi", resolution=0.5)

    # Get predicted labels
    predictions = scanvi_model.predict()
    adata.obs["scanvi_pred"] = predictions
    pred_counts = adata.obs["scanvi_pred"].value_counts()
    print(f"  scANVI predictions: {pred_counts.to_dict()}")

    print(f"  scANVI: done ({time.time()-t0:.0f}s)")
    scanvi_success = True
except Exception as e:
    print(f"  scANVI: FAILED ({e})")
    import traceback
    traceback.print_exc()
    scanvi_success = False

# ============================================================
# 3. Scanorama (fix)
# ============================================================
print("\n=== 3. Scanorama (fixed) ===")
try:
    import scanorama

    t0 = time.time()
    platforms = sorted(adata.obs["platform"].unique())
    adatas_by_platform = []
    for p in platforms:
        a = adata[adata.obs["platform"] == p].copy()
        a.X = a.layers["counts"].copy()
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        adatas_by_platform.append(a)

    # integrate_scanpy modifies adatas in-place and returns None
    scanorama.integrate_scanpy(adatas_by_platform, dimred=50)

    # Reassemble embeddings in original cell order
    emb_list = []
    cell_order = []
    for a in adatas_by_platform:
        if "X_scanorama" in a.obsm:
            emb_list.append(a.obsm["X_scanorama"])
            cell_order.extend(a.obs_names.tolist())
        else:
            print(f"  WARNING: X_scanorama not found in integrated adata (keys: {list(a.obsm.keys())})")

    if emb_list:
        emb_all = np.vstack(emb_list)
        # Reorder to match original adata
        cell_to_idx = {c: i for i, c in enumerate(cell_order)}
        reorder = [cell_to_idx[c] for c in adata.obs_names]
        adata.obsm["X_scanorama"] = emb_all[reorder].astype(np.float32)

        sc.pp.neighbors(adata, use_rep="X_scanorama", key_added="neighbors_scanorama")
        sc.tl.umap(adata, neighbors_key="neighbors_scanorama")
        adata.obsm["X_umap_scanorama"] = adata.obsm["X_umap"].copy()
        sc.tl.leiden(adata, neighbors_key="neighbors_scanorama", key_added="leiden_scanorama", resolution=0.5)
        print(f"  Scanorama: done ({time.time()-t0:.0f}s), embedding shape={adata.obsm['X_scanorama'].shape}")
        scanorama_success = True
    else:
        print("  Scanorama: no embeddings produced")
        scanorama_success = False
except Exception as e:
    print(f"  Scanorama: FAILED ({e})")
    import traceback
    traceback.print_exc()
    scanorama_success = False

# ============================================================
# 4. Per-cluster platform entropy
# ============================================================
print("\n=== 4. Per-cluster Platform Entropy ===")
from scipy.stats import entropy as scipy_entropy


def compute_platform_entropy(adata, cluster_key, platform_key="platform"):
    """Compute normalized Shannon entropy of platform distribution per cluster.
    1.0 = perfect mixing, 0.0 = single platform.
    """
    n_platforms = adata.obs[platform_key].nunique()
    max_ent = np.log(n_platforms)
    entropies = []
    for cluster in sorted(adata.obs[cluster_key].unique()):
        mask = adata.obs[cluster_key] == cluster
        counts = adata.obs.loc[mask, platform_key].value_counts()
        probs = counts / counts.sum()
        ent = scipy_entropy(probs) / max_ent if max_ent > 0 else 0
        entropies.append({"cluster": cluster, "entropy": ent, "n_cells": mask.sum()})
    return pd.DataFrame(entropies)


entropy_results = {}
for method_key in ["leiden_pca", "leiden_harmony", "leiden_scvi", "leiden_bbknn"]:
    if method_key not in adata.obs:
        continue
    method_name = method_key.replace("leiden_", "").upper()
    if method_name == "SCVI":
        method_name = "scVI"
    ent_df = compute_platform_entropy(adata, method_key)
    median_ent = ent_df["entropy"].median()
    entropy_results[method_name] = {"median_entropy": median_ent, "n_clusters": len(ent_df)}
    print(f"  {method_name}: median platform entropy = {median_ent:.3f} ({len(ent_df)} clusters)")
    print(f"    Criterion (>0.5): {'PASS' if median_ent > 0.5 else 'FAIL'}")

if scanvi_success and "leiden_scanvi" in adata.obs:
    ent_df = compute_platform_entropy(adata, "leiden_scanvi")
    median_ent = ent_df["entropy"].median()
    entropy_results["scANVI"] = {"median_entropy": median_ent, "n_clusters": len(ent_df)}
    print(f"  scANVI: median platform entropy = {median_ent:.3f} ({len(ent_df)} clusters)")
    print(f"    Criterion (>0.5): {'PASS' if median_ent > 0.5 else 'FAIL'}")

if scanorama_success and "leiden_scanorama" in adata.obs:
    ent_df = compute_platform_entropy(adata, "leiden_scanorama")
    median_ent = ent_df["entropy"].median()
    entropy_results["Scanorama"] = {"median_entropy": median_ent, "n_clusters": len(ent_df)}
    print(f"  Scanorama: median platform entropy = {median_ent:.3f} ({len(ent_df)} clusters)")

# ============================================================
# 5. Benchmark new methods
# ============================================================
print("\n=== 5. Benchmark (all methods including new) ===")
from sklearn.metrics import silhouette_score


def compute_batch_asw(adata, emb_key, batch_key="platform"):
    emb = adata.obsm[emb_key]
    labels = adata.obs[batch_key].values
    if len(np.unique(labels)) < 2:
        return np.nan
    return silhouette_score(emb, labels, sample_size=min(10000, adata.n_obs), random_state=42)


def compute_celltype_asw(adata, emb_key, ct_key="celltype"):
    if ct_key not in adata.obs:
        return np.nan
    labels = adata.obs[ct_key].values
    valid = ~pd.isna(labels)
    if valid.sum() < 100:
        return np.nan
    unique = np.unique(labels[valid])
    if len(unique) < 2:
        return np.nan
    return silhouette_score(adata[valid].obsm[emb_key], labels[valid],
                           sample_size=min(10000, valid.sum()), random_state=42)


all_methods = {
    "PCA": "X_pca",
    "Harmony": "X_harmony",
    "scVI": "X_scvi",
    "BBKNN": "X_pca",  # graph-based
}
if scanvi_success:
    all_methods["scANVI"] = "X_scanvi"
if scanorama_success:
    all_methods["Scanorama"] = "X_scanorama"

rows = []
for method, emb_key in all_methods.items():
    batch_asw = compute_batch_asw(adata, emb_key)
    ct_asw = compute_celltype_asw(adata, emb_key, "celltype")
    ct_broad_asw = compute_celltype_asw(adata, emb_key, "celltype_broad")
    ent_info = entropy_results.get(method, {})

    row = {
        "method": method,
        "batch_asw": batch_asw,
        "celltype_asw": ct_asw,
        "celltype_broad_asw": ct_broad_asw,
        "batch_score": 1 - abs(batch_asw),
        "platform_entropy": ent_info.get("median_entropy", np.nan),
    }
    rows.append(row)
    print(f"  {method}: batch_ASW={batch_asw:.3f}, ct_broad_ASW={ct_broad_asw:.3f}, "
          f"batch_score={row['batch_score']:.3f}, entropy={row['platform_entropy']:.3f}")

benchmark_df = pd.DataFrame(rows).sort_values("batch_score", ascending=False)
benchmark_df.to_csv(OUTDIR / "m2_followup_benchmark.csv", index=False)
print(f"\nBenchmark saved to {OUTDIR / 'm2_followup_benchmark.csv'}")
print(benchmark_df.to_string(index=False))

# ============================================================
# 6. scANVI-specific evaluation
# ============================================================
if scanvi_success:
    print("\n=== 6. scANVI-specific evaluation ===")

    # Label transfer accuracy: compare scANVI predictions to original labels
    ct_key = "celltype_broad"
    orig = adata.obs[ct_key].values
    pred = adata.obs["scanvi_pred"].values
    valid = ~pd.isna(orig) & (pred != "Unknown")
    if valid.sum() > 0:
        acc = (orig[valid] == pred[valid]).mean()
        print(f"  scANVI prediction accuracy vs original {ct_key}: {acc:.3f} ({valid.sum()} cells)")

        # Per-platform accuracy
        for platform in sorted(adata.obs["platform"].unique()):
            pmask = (adata.obs["platform"] == platform).values & valid
            if pmask.sum() > 0:
                pacc = (orig[pmask] == pred[pmask]).mean()
                print(f"    {platform}: {pacc:.3f} ({pmask.sum()} cells)")

    # Cross-platform kNN transfer with scANVI embedding
    from sklearn.neighbors import KNeighborsClassifier

    print("\n  Cross-platform kNN transfer (scANVI embedding):")
    for ct_key in ["celltype", "celltype_broad"]:
        labels = adata.obs[ct_key].values
        valid = ~pd.isna(labels)
        cosmx_mask = (adata.obs["platform"] == "CosMx").values & valid
        xenium_mask = (adata.obs["platform"] == "Xenium").values & valid
        if cosmx_mask.sum() < 20 or xenium_mask.sum() < 20:
            continue
        emb = adata.obsm["X_scanvi"]
        knn = KNeighborsClassifier(n_neighbors=15, metric="euclidean")
        knn.fit(emb[cosmx_mask], labels[cosmx_mask])
        acc_c2x = knn.score(emb[xenium_mask], labels[xenium_mask])
        knn2 = KNeighborsClassifier(n_neighbors=15, metric="euclidean")
        knn2.fit(emb[xenium_mask], labels[xenium_mask])
        acc_x2c = knn2.score(emb[cosmx_mask], labels[cosmx_mask])
        print(f"    {ct_key}: CosMx->Xenium={acc_c2x:.3f}, Xenium->CosMx={acc_x2c:.3f}")

# ============================================================
# 7. Formal success criteria (from plan)
# ============================================================
print("\n=== 7. Formal Success Criteria (from plan) ===")
print("| Criterion | Threshold | M2 Result | Status |")
print("|-----------|-----------|-----------|--------|")

# ASW batch
best = benchmark_df.iloc[0]
asw_pass = abs(best["batch_asw"]) < 0.5
print(f"| ASW_batch (best method) | < 0.5 | {best['batch_asw']:.3f} ({best['method']}) | {'PASS' if asw_pass else 'FAIL'} |")

# Platform entropy
best_ent = max(entropy_results.values(), key=lambda x: x["median_entropy"]) if entropy_results else {}
ent_val = best_ent.get("median_entropy", 0)
ent_pass = ent_val > 0.5
ent_method = [k for k, v in entropy_results.items() if v.get("median_entropy") == ent_val][0] if entropy_results else "N/A"
print(f"| Per-cluster platform entropy | Median > 0.5 | {ent_val:.3f} ({ent_method}) | {'PASS' if ent_pass else 'FAIL'} |")

# Bio conservation
bio_pass = best["celltype_broad_asw"] >= 0
print(f"| Bio conservation (ct_broad ASW >= 0) | >= 0 | {best['celltype_broad_asw']:.3f} ({best['method']}) | {'PASS' if bio_pass else 'FAIL'} |")

# IBD markers
marker_pass = len(present) >= 5
print(f"| IBD markers in intersection | >= 5/8 | {len(present)}/8 | {'PASS' if marker_pass else 'FAIL'} |")

# Diagnosis signal (UC vs Healthy — no CD in M2)
print(f"| Diagnosis signal (CD vs UC) | Fisher p < 0.05 | N/A (no CD in M2) | N/A |")

# ============================================================
# 8. UMAP plots for new methods
# ============================================================
print("\n=== 8. Plotting ===")
color_keys = ["platform", "patient_id"]
if "disease" in adata.obs:
    color_keys.append("disease")
if "celltype_broad" in adata.obs:
    color_keys.append("celltype_broad")
if scanvi_success and "scanvi_pred" in adata.obs:
    color_keys.append("scanvi_pred")

n_colors = len(color_keys)

new_methods = {}
if scanvi_success:
    new_methods["scANVI"] = "X_umap_scanvi"
if scanorama_success:
    new_methods["Scanorama"] = "X_umap_scanorama"

for method, umap_key in new_methods.items():
    if umap_key not in adata.obsm:
        continue
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

# Updated UMAP grid with all methods
all_method_umaps = ["PCA", "Harmony", "scVI", "BBKNN"]
if scanvi_success:
    all_method_umaps.append("scANVI")
if scanorama_success:
    all_method_umaps.append("Scanorama")

n_methods = len(all_method_umaps)
fig, axes = plt.subplots(n_methods, n_colors, figsize=(6 * n_colors, 5 * n_methods))
if n_methods == 1:
    axes = axes.reshape(1, -1)

for i, method in enumerate(all_method_umaps):
    umap_key = f"X_umap_{method.lower()}"
    if umap_key not in adata.obsm:
        umap_key = "X_umap_pca"
    adata.obsm["X_umap"] = adata.obsm[umap_key]
    for j, key in enumerate(color_keys):
        sc.pl.umap(adata, color=key, ax=axes[i, j], show=False, title=f"{method} - {key}")

plt.tight_layout()
fig.savefig(FIGDIR / "m2_umap_grid_all.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved: m2_umap_grid_all.png")

# ============================================================
# 9. Save updated adata
# ============================================================
print("\n=== Saving ===")
adata.write_h5ad(OUTDIR / "adata.m2.h5ad")
print(f"Saved: {OUTDIR / 'adata.m2.h5ad'} ({adata.n_obs} cells x {adata.n_vars} genes)")

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 60)
print("M2 FOLLOW-UP COMPLETE")
print("=" * 60)
print(f"IBD markers: {len(present)}/8 present")
print(f"Methods run: {list(all_methods.keys())}")
print(f"Best method: {best['method']} (batch_score={best['batch_score']:.3f})")
if scanvi_success:
    scanvi_row = benchmark_df[benchmark_df["method"] == "scANVI"].iloc[0]
    print(f"scANVI: batch_score={scanvi_row['batch_score']:.3f}, ct_broad_ASW={scanvi_row['celltype_broad_asw']:.3f}")
