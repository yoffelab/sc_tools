#!/usr/bin/env python3
"""M1 follow-up: proper biological evaluation of integration quality.

Addresses issues with naive celltype ASW:
1. Per-platform celltype ASW (are labels self-consistent in 119-gene space?)
2. celltype_broad as primary bio metric
3. Re-cluster on shared genes, compare to original labels via ARI/NMI
4. Cross-platform label transfer accuracy (kNN classifier)

Reads: results/m1_benchmark/adata.m1.h5ad
Writes: results/m1_benchmark/m1_bio_eval.csv, figures/m1/m1_bio_*.png
"""

import sys
import os
from pathlib import Path

# Force unbuffered output
os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout.reconfigure(line_buffering=True)

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score, adjusted_rand_score, normalized_mutual_info_score
from sklearn.neighbors import KNeighborsClassifier

WORKDIR = Path("/home/fs01/juk4007/elementolab/projects/ibd_spatial")
OUTDIR = WORKDIR / "results" / "m1_benchmark"
FIGDIR = WORKDIR / "figures" / "m1"

MAX_CELLS = 50000  # subsample for faster metrics

print("=== Loading M1 adata ===")
adata = ad.read_h5ad(OUTDIR / "adata.m1.h5ad")
print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

# Subsample for metrics (stratified by platform)
if adata.n_obs > MAX_CELLS:
    np.random.seed(42)
    idx = []
    for platform in adata.obs["platform"].unique():
        pmask = adata.obs["platform"] == platform
        pidx = np.where(pmask)[0]
        n_take = min(len(pidx), MAX_CELLS // 2)
        idx.extend(np.random.choice(pidx, n_take, replace=False))
    idx = sorted(idx)
    adata = adata[idx].copy()
    print(f"Subsampled to {adata.n_obs} cells for metrics")

# Identify integration embeddings
emb_keys = {
    "PCA": "X_pca",
    "Harmony": "X_harmony",
    "scVI": "X_scvi",
}
# Only include methods that exist
emb_keys = {k: v for k, v in emb_keys.items() if v in adata.obsm}
print(f"Methods: {list(emb_keys.keys())}")

# ============================================================
# 1. Per-platform celltype ASW (label self-consistency check)
# ============================================================
print("\n=== 1. Per-platform celltype ASW ===")
print("Do celltype labels form coherent clusters in the 119-gene embedding?")

ct_keys = []
if "celltype" in adata.obs.columns:
    ct_keys.append("celltype")
if "celltype_broad" in adata.obs.columns:
    ct_keys.append("celltype_broad")

platform_asw_rows = []
for platform in sorted(adata.obs["platform"].unique()):
    mask = adata.obs["platform"] == platform
    sub = adata[mask]
    for ct_key in ct_keys:
        labels = sub.obs[ct_key].values
        valid = ~pd.isna(labels)
        if valid.sum() < 100:
            continue
        unique_labels = np.unique(labels[valid])
        if len(unique_labels) < 2:
            continue
        for method, emb_key in emb_keys.items():
            emb = sub[valid].obsm[emb_key]
            asw = silhouette_score(
                emb, labels[valid],
                sample_size=min(5000, valid.sum()),
                random_state=42,
            )
            platform_asw_rows.append({
                "platform": platform,
                "celltype_key": ct_key,
                "method": method,
                "asw": asw,
                "n_cells": valid.sum(),
                "n_types": len(unique_labels),
            })
            print(f"  {platform} / {ct_key} / {method}: ASW={asw:.3f} ({len(unique_labels)} types, {valid.sum()} cells)")

# Cross-platform (original metric for comparison)
for ct_key in ct_keys:
    labels = adata.obs[ct_key].values
    valid = ~pd.isna(labels)
    if valid.sum() < 100:
        continue
    unique_labels = np.unique(labels[valid])
    if len(unique_labels) < 2:
        continue
    for method, emb_key in emb_keys.items():
        emb = adata[valid].obsm[emb_key]
        asw = silhouette_score(
            emb, labels[valid],
            sample_size=min(5000, valid.sum()),
            random_state=42,
        )
        platform_asw_rows.append({
            "platform": "Combined",
            "celltype_key": ct_key,
            "method": method,
            "asw": asw,
            "n_cells": valid.sum(),
            "n_types": len(unique_labels),
        })
        print(f"  Combined / {ct_key} / {method}: ASW={asw:.3f}")

platform_asw_df = pd.DataFrame(platform_asw_rows)

# ============================================================
# 2. Re-cluster and compare to original labels (ARI/NMI)
# ============================================================
print("\n=== 2. Re-clustering vs original labels (ARI/NMI) ===")

cluster_rows = []
for method, emb_key in emb_keys.items():
    neighbors_key = f"neighbors_{method.lower()}"
    if neighbors_key not in adata.uns:
        # Build neighbors
        sc.pp.neighbors(adata, use_rep=emb_key, key_added=neighbors_key)

    for res in [0.3, 0.5, 1.0]:
        leiden_key = f"leiden_{method.lower()}_r{res}"
        sc.tl.leiden(
            adata,
            neighbors_key=neighbors_key,
            key_added=leiden_key,
            resolution=res,
        )
        n_clusters = adata.obs[leiden_key].nunique()

        for ct_key in ct_keys:
            labels = adata.obs[ct_key].values
            valid = ~pd.isna(labels)
            if valid.sum() < 100:
                continue
            clusters = adata.obs[leiden_key].values[valid]
            ari = adjusted_rand_score(labels[valid], clusters)
            nmi = normalized_mutual_info_score(labels[valid], clusters)
            cluster_rows.append({
                "method": method,
                "resolution": res,
                "n_clusters": n_clusters,
                "celltype_key": ct_key,
                "ARI": ari,
                "NMI": nmi,
            })
            print(f"  {method} res={res} ({n_clusters} clusters) vs {ct_key}: ARI={ari:.3f}, NMI={nmi:.3f}")

cluster_df = pd.DataFrame(cluster_rows)

# ============================================================
# 3. Cross-platform label transfer (kNN classifier)
# ============================================================
print("\n=== 3. Cross-platform label transfer (kNN) ===")
print("Train on one platform, predict on other — tests if integration aligns cell types")

transfer_rows = []
for ct_key in ct_keys:
    labels = adata.obs[ct_key].values
    valid = ~pd.isna(labels)
    if valid.sum() < 100:
        continue

    cosmx_mask = (adata.obs["platform"] == "CosMx").values & valid
    xenium_mask = (adata.obs["platform"] == "Xenium").values & valid

    if cosmx_mask.sum() < 50 or xenium_mask.sum() < 50:
        continue

    for method, emb_key in emb_keys.items():
        emb = adata.obsm[emb_key]

        # Train on CosMx, predict Xenium
        knn = KNeighborsClassifier(n_neighbors=15, metric="euclidean")
        knn.fit(emb[cosmx_mask], labels[cosmx_mask])
        acc_c2x = knn.score(emb[xenium_mask], labels[xenium_mask])

        # Train on Xenium, predict CosMx
        knn2 = KNeighborsClassifier(n_neighbors=15, metric="euclidean")
        knn2.fit(emb[xenium_mask], labels[xenium_mask])
        acc_x2c = knn2.score(emb[cosmx_mask], labels[cosmx_mask])

        transfer_rows.append({
            "method": method,
            "celltype_key": ct_key,
            "acc_cosmx_to_xenium": acc_c2x,
            "acc_xenium_to_cosmx": acc_x2c,
            "mean_transfer_acc": (acc_c2x + acc_x2c) / 2,
        })
        print(f"  {method} / {ct_key}: CosMx→Xenium={acc_c2x:.3f}, Xenium→CosMx={acc_x2c:.3f}, mean={((acc_c2x+acc_x2c)/2):.3f}")

transfer_df = pd.DataFrame(transfer_rows)

# ============================================================
# 4. Batch ASW within cell types (mixing per cell type)
# ============================================================
print("\n=== 4. Per-celltype platform mixing ===")
print("Within each cell type, do platforms mix? (ASW of platform label)")

perct_rows = []
ct_key = "celltype_broad" if "celltype_broad" in adata.obs.columns else "celltype"
labels = adata.obs[ct_key].values
valid = ~pd.isna(labels)

for ct in sorted(adata.obs.loc[valid, ct_key].unique()):
    ct_mask = valid & (adata.obs[ct_key] == ct).values
    if ct_mask.sum() < 50:
        continue
    sub = adata[ct_mask]
    if sub.obs["platform"].nunique() < 2:
        continue

    for method, emb_key in emb_keys.items():
        emb = sub.obsm[emb_key]
        asw = silhouette_score(
            emb, sub.obs["platform"].values,
            sample_size=min(2000, sub.n_obs),
            random_state=42,
        )
        perct_rows.append({
            "celltype": ct,
            "method": method,
            "platform_asw": asw,
            "n_cells": sub.n_obs,
        })
    print(f"  {ct}: n={ct_mask.sum()}")

perct_df = pd.DataFrame(perct_rows)
if len(perct_df) > 0:
    # Pivot for display
    pivot = perct_df.pivot_table(index="celltype", columns="method", values="platform_asw")
    print(pivot.round(3).to_string())

# ============================================================
# 5. Summary plots
# ============================================================
print("\n=== Plotting ===")

# 5a. Bar plot: per-platform celltype ASW
if len(platform_asw_df) > 0:
    fig, axes = plt.subplots(1, len(ct_keys), figsize=(7 * len(ct_keys), 5))
    if len(ct_keys) == 1:
        axes = [axes]
    for i, ct_key in enumerate(ct_keys):
        sub = platform_asw_df[platform_asw_df["celltype_key"] == ct_key]
        if len(sub) == 0:
            continue
        pivot = sub.pivot_table(index="platform", columns="method", values="asw")
        pivot.plot(kind="bar", ax=axes[i], rot=0)
        axes[i].set_title(f"Celltype ASW ({ct_key})")
        axes[i].set_ylabel("Silhouette Score")
        axes[i].axhline(0, color="k", linestyle="--", alpha=0.3)
    plt.tight_layout()
    fig.savefig(FIGDIR / "m1_bio_platform_asw.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved: m1_bio_platform_asw.png")

# 5b. ARI/NMI bar plot
if len(cluster_df) > 0:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for i, metric in enumerate(["ARI", "NMI"]):
        sub = cluster_df[cluster_df["celltype_key"] == ct_keys[0]]
        pivot = sub.pivot_table(index="resolution", columns="method", values=metric)
        pivot.plot(kind="bar", ax=axes[i], rot=0)
        axes[i].set_title(f"{metric} (leiden vs {ct_keys[0]})")
        axes[i].set_ylabel(metric)
    plt.tight_layout()
    fig.savefig(FIGDIR / "m1_bio_ari_nmi.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved: m1_bio_ari_nmi.png")

# 5c. Transfer accuracy bar plot
if len(transfer_df) > 0:
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(transfer_df["method"].unique()))
    width = 0.25
    for i, ct_key in enumerate(ct_keys):
        sub = transfer_df[transfer_df["celltype_key"] == ct_key]
        ax.bar(x + i * width, sub["mean_transfer_acc"].values, width, label=ct_key)
    ax.set_xticks(x + width * (len(ct_keys) - 1) / 2)
    ax.set_xticklabels(transfer_df["method"].unique())
    ax.set_ylabel("Mean transfer accuracy")
    ax.set_title("Cross-platform label transfer (kNN, k=15)")
    ax.legend()
    ax.set_ylim(0, 1)
    plt.tight_layout()
    fig.savefig(FIGDIR / "m1_bio_transfer.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved: m1_bio_transfer.png")

# ============================================================
# 6. Save all results
# ============================================================
print("\n=== Saving ===")
platform_asw_df.to_csv(OUTDIR / "m1_bio_platform_asw.csv", index=False)
cluster_df.to_csv(OUTDIR / "m1_bio_cluster_vs_labels.csv", index=False)
if len(transfer_df) > 0:
    transfer_df.to_csv(OUTDIR / "m1_bio_label_transfer.csv", index=False)
if len(perct_df) > 0:
    perct_df.to_csv(OUTDIR / "m1_bio_perct_mixing.csv", index=False)
print(f"Saved CSVs to {OUTDIR}")

# ============================================================
# 7. Summary
# ============================================================
print("\n" + "=" * 60)
print("M1 BIO EVALUATION SUMMARY")
print("=" * 60)

print("\n1. Per-platform celltype ASW (label self-consistency in 119-gene space):")
for platform in ["CosMx", "Xenium", "Combined"]:
    sub = platform_asw_df[(platform_asw_df["platform"] == platform) & (platform_asw_df["celltype_key"] == ct_keys[0])]
    if len(sub) > 0:
        best = sub.loc[sub["asw"].idxmax()]
        print(f"   {platform}: best={best['method']} ASW={best['asw']:.3f}")

print(f"\n2. Best ARI/NMI (leiden re-clustering vs {ct_keys[0]}):")
if len(cluster_df) > 0:
    best_ari = cluster_df.loc[cluster_df["ARI"].idxmax()]
    best_nmi = cluster_df.loc[cluster_df["NMI"].idxmax()]
    print(f"   Best ARI: {best_ari['method']} res={best_ari['resolution']} ARI={best_ari['ARI']:.3f}")
    print(f"   Best NMI: {best_nmi['method']} res={best_nmi['resolution']} NMI={best_nmi['NMI']:.3f}")

print("\n3. Cross-platform label transfer accuracy:")
if len(transfer_df) > 0:
    for _, row in transfer_df.iterrows():
        print(f"   {row['method']} / {row['celltype_key']}: {row['mean_transfer_acc']:.3f}")

print("\n4. Interpretation:")
# Check if per-platform ASW is also negative
cosmx_asw = platform_asw_df[
    (platform_asw_df["platform"] == "CosMx") &
    (platform_asw_df["celltype_key"] == ct_keys[0])
]["asw"].values
if len(cosmx_asw) > 0 and np.mean(cosmx_asw) < 0:
    print("   -> Per-platform celltype ASW also negative: 119 genes cannot resolve fine cell types")
    print("      even WITHIN a single platform. The negative combined ASW is NOT an integration artifact.")
elif len(cosmx_asw) > 0:
    print("   -> Per-platform celltype ASW is positive: labels are self-consistent within platforms.")
    print("      Negative combined ASW suggests integration is mixing cell types across platforms.")
