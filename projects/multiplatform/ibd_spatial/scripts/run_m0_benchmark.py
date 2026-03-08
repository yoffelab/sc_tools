#!/usr/bin/env python3
"""Milestone 0: noseg vs withseg integration benchmark.

Scope: 4 noseg + 4 withseg samples (same 4 patients, Rectum, 377-gene Xenium panel).
This is the technical replicate baseline — expect near-perfect integration.

Usage:
    python run_m0_benchmark.py

Reads from: data/xenium_noseg_*/adata.p0.h5ad + data/xenium_withseg_*/adata.p0.h5ad
Writes to:  results/m0_benchmark/
"""

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

WORKDIR = Path("/home/fs01/juk4007/elementolab/projects/ibd_spatial")
OUTDIR = WORKDIR / "results" / "m0_benchmark"
FIGDIR = WORKDIR / "figures" / "m0"

OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)

# --- 1. Load and concatenate ---
print("=== Loading samples ===")
adatas = []
for panel in ["xenium_noseg", "xenium_withseg"]:
    for i in range(1, 5):
        sample_id = f"{panel}_{i:02d}"
        path = WORKDIR / "data" / sample_id / "adata.p0.h5ad"
        if not path.exists():
            print(f"  MISSING: {path}")
            continue
        a = ad.read_h5ad(path)
        # Ensure consistent obs columns
        a.obs["panel_variant"] = panel  # noseg vs withseg
        a.obs["sample"] = sample_id
        a.var_names_make_unique()
        adatas.append(a)
        print(f"  {sample_id}: {a.shape[0]} cells x {a.shape[1]} genes")

adata = ad.concat(adatas, join="inner", merge="same")
adata.obs_names_make_unique()
print(f"\nConcatenated: {adata.shape[0]} cells x {adata.shape[1]} genes")
print(f"Panels: {adata.obs['panel_variant'].value_counts().to_dict()}")
print(f"Patients: {adata.obs['patient_id'].value_counts().to_dict()}")

# --- 2. Basic QC ---
print("\n=== QC ===")
sc.pp.calculate_qc_metrics(adata, percent_top=[50, 100, 200], inplace=True)
print(f"Median counts/cell: {adata.obs['total_counts'].median():.0f}")
print(f"Median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}")

# Filter low-quality cells
min_counts = 10
min_genes = 5
n_before = adata.n_obs
sc.pp.filter_cells(adata, min_counts=min_counts)
sc.pp.filter_cells(adata, min_genes=min_genes)
print(f"After filtering (min_counts={min_counts}, min_genes={min_genes}): {adata.n_obs} cells (removed {n_before - adata.n_obs})")

# Filter genes expressed in < 10 cells
sc.pp.filter_genes(adata, min_cells=10)
print(f"Genes after filtering: {adata.n_vars}")

# --- 3. Preprocessing ---
print("\n=== Preprocessing ===")
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# For targeted panels: use all genes (no HVG selection)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=min(50, adata.n_vars - 1))
print(f"PCA: {adata.obsm['X_pca'].shape}")

# Unintegrated baseline
sc.pp.neighbors(adata, use_rep="X_pca", key_added="neighbors_pca")
sc.tl.umap(adata, neighbors_key="neighbors_pca")
adata.obsm["X_umap_pca"] = adata.obsm["X_umap"].copy()
sc.tl.leiden(adata, neighbors_key="neighbors_pca", key_added="leiden_pca", resolution=0.5)

# --- 4. Integration methods ---
print("\n=== Integration ===")
results = {}

# 4a. Unintegrated PCA (baseline)
results["PCA"] = {"embedding_key": "X_pca", "neighbors_key": "neighbors_pca"}

# 4b. Harmony
try:
    import harmonypy

    adata_harmony = adata.copy()
    ho = harmonypy.run_harmony(adata_harmony.obsm["X_pca"], adata_harmony.obs, "panel_variant")
    emb = ho.Z_corr
    # Handle PyTorch tensor
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
    print("  Harmony: done")
except Exception as e:
    print(f"  Harmony: FAILED ({e})")

# 4c. scVI (requires raw counts)
try:
    import scvi

    adata_scvi = adata.copy()
    adata_scvi.X = adata_scvi.layers["counts"].copy()
    scvi.model.SCVI.setup_anndata(adata_scvi, batch_key="panel_variant")
    model = scvi.model.SCVI(adata_scvi, n_latent=10, n_hidden=64, n_layers=1)
    model.train(max_epochs=100, early_stopping=True, batch_size=256)
    latent = model.get_latent_representation()
    adata.obsm["X_scvi"] = latent.astype(np.float32)
    sc.pp.neighbors(adata, use_rep="X_scvi", key_added="neighbors_scvi")
    sc.tl.umap(adata, neighbors_key="neighbors_scvi")
    adata.obsm["X_umap_scvi"] = adata.obsm["X_umap"].copy()
    sc.tl.leiden(adata, neighbors_key="neighbors_scvi", key_added="leiden_scvi", resolution=0.5)
    results["scVI"] = {"embedding_key": "X_scvi", "neighbors_key": "neighbors_scvi"}
    print("  scVI: done")
except Exception as e:
    print(f"  scVI: FAILED ({e})")

# 4d. BBKNN
try:
    import bbknn

    adata_bbknn = adata.copy()
    bbknn.bbknn(adata_bbknn, batch_key="panel_variant", use_rep="X_pca")
    # BBKNN modifies connectivities in place
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
    print("  BBKNN: done")
except Exception as e:
    print(f"  BBKNN: FAILED ({e})")

# 4e. Scanorama
try:
    import scanorama

    panels = adata.obs["panel_variant"].cat.categories if hasattr(adata.obs["panel_variant"], "cat") else adata.obs["panel_variant"].unique()
    adatas_by_panel = [adata[adata.obs["panel_variant"] == p].copy() for p in panels]
    corrected, _ = scanorama.correct_scanpy(adatas_by_panel)
    # Reassemble
    corrected_X = np.vstack([a.obsm["X_scanorama"] for a in corrected])
    adata.obsm["X_scanorama"] = corrected_X.astype(np.float32)
    sc.pp.neighbors(adata, use_rep="X_scanorama", key_added="neighbors_scanorama")
    sc.tl.umap(adata, neighbors_key="neighbors_scanorama")
    adata.obsm["X_umap_scanorama"] = adata.obsm["X_umap"].copy()
    sc.tl.leiden(adata, neighbors_key="neighbors_scanorama", key_added="leiden_scanorama", resolution=0.5)
    results["Scanorama"] = {"embedding_key": "X_scanorama", "neighbors_key": "neighbors_scanorama"}
    print("  Scanorama: done")
except Exception as e:
    print(f"  Scanorama: FAILED ({e})")

# --- 5. Benchmark metrics ---
print("\n=== Benchmark ===")
from sklearn.metrics import silhouette_score


def compute_batch_asw(adata, embedding_key, batch_key="panel_variant"):
    """ASW for batch: lower = better mixing (negate silhouette)."""
    emb = adata.obsm[embedding_key]
    labels = adata.obs[batch_key].values
    if len(np.unique(labels)) < 2:
        return np.nan
    asw = silhouette_score(emb, labels, sample_size=min(5000, adata.n_obs), random_state=42)
    return asw


def compute_patient_mixing(adata, embedding_key, patient_key="patient_id", batch_key="panel_variant"):
    """Per-patient cross-platform mixing score.
    For each patient, compute silhouette of batch labels within that patient's cells.
    Perfect mixing = 0; perfect separation = 1.
    We want this close to 0 (platforms mixed within same patient).
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


benchmark_rows = []
for method, info in results.items():
    emb_key = info["embedding_key"]
    batch_asw = compute_batch_asw(adata, emb_key)
    patient_mix = compute_patient_mixing(adata, emb_key)

    row = {
        "method": method,
        "batch_asw": batch_asw,
        "patient_mixing_asw": patient_mix,
        "batch_score": 1 - abs(batch_asw),  # higher = better mixing
    }
    benchmark_rows.append(row)
    print(f"  {method}: batch_ASW={batch_asw:.3f}, patient_mix={patient_mix:.3f}, batch_score={row['batch_score']:.3f}")

benchmark_df = pd.DataFrame(benchmark_rows).sort_values("batch_score", ascending=False)
benchmark_df.to_csv(OUTDIR / "m0_benchmark.csv", index=False)
print(f"\nBenchmark saved to {OUTDIR / 'm0_benchmark.csv'}")
print(benchmark_df.to_string(index=False))

# --- 6. UMAP plots ---
print("\n=== Plotting ===")
for method, info in results.items():
    umap_key = f"X_umap_{method.lower()}"
    if umap_key not in adata.obsm:
        umap_key = "X_umap_pca"

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f"M0: {method}", fontsize=14, fontweight="bold")

    # Color by panel
    adata.obsm["X_umap"] = adata.obsm[umap_key]
    sc.pl.umap(adata, color="panel_variant", ax=axes[0], show=False, title="Panel variant")
    sc.pl.umap(adata, color="patient_id", ax=axes[1], show=False, title="Patient")
    sc.pl.umap(adata, color="disease", ax=axes[2], show=False, title="Disease")

    plt.tight_layout()
    fig.savefig(FIGDIR / f"m0_umap_{method.lower()}.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: m0_umap_{method.lower()}.png")

# Summary UMAP grid
n_methods = len(results)
fig, axes = plt.subplots(n_methods, 3, figsize=(18, 5 * n_methods))
if n_methods == 1:
    axes = axes.reshape(1, -1)

for i, (method, info) in enumerate(results.items()):
    umap_key = f"X_umap_{method.lower()}"
    if umap_key not in adata.obsm:
        umap_key = "X_umap_pca"
    adata.obsm["X_umap"] = adata.obsm[umap_key]

    sc.pl.umap(adata, color="panel_variant", ax=axes[i, 0], show=False, title=f"{method} — Panel")
    sc.pl.umap(adata, color="patient_id", ax=axes[i, 1], show=False, title=f"{method} — Patient")
    sc.pl.umap(adata, color="disease", ax=axes[i, 2], show=False, title=f"{method} — Disease")

plt.tight_layout()
fig.savefig(FIGDIR / "m0_umap_grid.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved: m0_umap_grid.png")

# Save concatenated adata
adata.write_h5ad(OUTDIR / "adata.m0.h5ad")
print(f"\nSaved: {OUTDIR / 'adata.m0.h5ad'} ({adata.n_obs} cells x {adata.n_vars} genes)")

# --- 7. Summary ---
print("\n" + "=" * 60)
print("MILESTONE 0 COMPLETE")
print("=" * 60)
best = benchmark_df.iloc[0]
print(f"Best method: {best['method']} (batch_score={best['batch_score']:.3f})")
print(f"Expected: near-perfect mixing (batch_score > 0.8) since same panel/patients")
if best["batch_score"] < 0.8:
    print("WARNING: batch_score < 0.8 — unexpected for technical replicates")
print(f"\nOutputs:")
print(f"  Benchmark: {OUTDIR / 'm0_benchmark.csv'}")
print(f"  Figures: {FIGDIR}")
print(f"  AnnData: {OUTDIR / 'adata.m0.h5ad'}")
