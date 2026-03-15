#!/usr/bin/env python3
"""Run resolVI with semisupervised=True for M1 and M2.

resolVI natively supports semi-supervised training via a Gaussian mixture
prior per cell type + classifier in latent space. This combines spatial
awareness (resolVI) with cell type label propagation (like scANVI) in a
single model.

Compares against:
  - Existing unsupervised resolVI (X_resolvi)
  - Existing scANVI (X_scanvi)

Reads:
  results/m1_benchmark/adata.m1.h5ad
  results/m2_benchmark/adata.m2.h5ad
Writes:
  Updated adata with X_resolvi_ss embeddings
  results/resolvi_semisup_benchmark.csv
  figures/m{1,2}/m{1,2}_umap_resolvi_ss.png
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

import torch
torch.set_float32_matmul_precision("medium")

from scvi.external import RESOLVI

WORKDIR = Path("/home/fs01/juk4007/elementolab/projects/ibd_spatial")


def run_resolvi_semisup(adata, batch_key, label_key, spatial_key="spatial",
                        n_latent=20, max_epochs=100):
    """Run resolVI with semisupervised=True."""
    print(f"  Setting up resolVI semi-supervised: batch_key={batch_key}, "
          f"label_key={label_key}, n_latent={n_latent}")

    adata_rv = adata.copy()

    # Ensure raw counts in X
    if "counts" in adata_rv.layers:
        adata_rv.X = adata_rv.layers["counts"].copy()

    # resolVI expects obsm['X_spatial']
    if "X_spatial" not in adata_rv.obsm and "spatial" in adata_rv.obsm:
        adata_rv.obsm["X_spatial"] = adata_rv.obsm["spatial"].copy()
        print("  Copied obsm['spatial'] -> obsm['X_spatial']")

    # HVG selection for high-gene datasets
    if adata_rv.n_vars > 500:
        sc.pp.highly_variable_genes(
            adata_rv, n_top_genes=2000, batch_key=batch_key, flavor="seurat_v3"
        )
        adata_rv = adata_rv[:, adata_rv.var["highly_variable"]].copy()
        print(f"  HVG subset: {adata_rv.n_vars} genes")

    # Prepare labels — handle NaN/missing
    adata_rv.obs["_resolvi_labels"] = (
        adata_rv.obs[label_key].astype(str).replace("nan", "Unknown")
    )
    n_labeled = (adata_rv.obs["_resolvi_labels"] != "Unknown").sum()
    n_types = adata_rv.obs["_resolvi_labels"].nunique()
    print(f"  Labels: {n_labeled}/{adata_rv.n_obs} cells labeled, {n_types} types")

    # Setup with labels_key for semi-supervised
    RESOLVI.setup_anndata(
        adata_rv,
        batch_key=batch_key,
        labels_key="_resolvi_labels",
        layer=None,
    )

    model = RESOLVI(
        adata_rv,
        n_latent=n_latent,
        semisupervised=True,
    )

    print(f"  Training resolVI semi-supervised ({max_epochs} epochs)...")
    t0 = time.time()
    model.train(
        max_epochs=max_epochs,
        batch_size=512,
    )
    elapsed = time.time() - t0
    print(f"  Training done ({elapsed:.0f}s)")

    return model, adata_rv


def compute_metrics(adata, emb_key, batch_key, ct_key=None):
    """Compute batch ASW and celltype ASW."""
    from sklearn.metrics import silhouette_score
    emb = adata.obsm[emb_key]
    n_sample = min(10000, adata.n_obs)

    batch_asw = silhouette_score(
        emb, adata.obs[batch_key].values,
        sample_size=n_sample, random_state=42
    )

    ct_asw = None
    ct_broad_asw = None
    if ct_key and ct_key in adata.obs:
        labels = adata.obs[ct_key].values
        valid = ~pd.isna(labels)
        if valid.sum() > 100 and len(np.unique(labels[valid])) > 1:
            ct_asw = silhouette_score(
                emb[valid], labels[valid],
                sample_size=min(n_sample, valid.sum()), random_state=42
            )

    ct_broad_key = "celltype_broad"
    if ct_broad_key in adata.obs:
        labels = adata.obs[ct_broad_key].values
        valid = ~pd.isna(labels)
        if valid.sum() > 100 and len(np.unique(labels[valid])) > 1:
            ct_broad_asw = silhouette_score(
                emb[valid], labels[valid],
                sample_size=min(n_sample, valid.sum()), random_state=42
            )

    return batch_asw, ct_asw, ct_broad_asw


def compute_platform_entropy(adata, emb_key, batch_key, resolution=0.5):
    """Compute per-cluster platform entropy."""
    from scipy.stats import entropy as scipy_entropy

    sc.pp.neighbors(adata, use_rep=emb_key, key_added="neighbors_tmp")
    sc.tl.leiden(adata, neighbors_key="neighbors_tmp", key_added="leiden_tmp",
                 resolution=resolution)

    n_batches = adata.obs[batch_key].nunique()
    max_ent = np.log(n_batches) if n_batches > 1 else 1.0
    entropies = []
    for cluster in adata.obs["leiden_tmp"].unique():
        mask = adata.obs["leiden_tmp"] == cluster
        counts = adata.obs.loc[mask, batch_key].value_counts()
        ent = scipy_entropy(counts / counts.sum()) / max_ent if max_ent > 0 else 0
        entropies.append(ent)
    return np.median(entropies), len(entropies)


def plot_umaps(adata, emb_key, method_name, milestone, figdir, batch_key, color_keys):
    """Generate UMAP from embedding and save."""
    sc.pp.neighbors(adata, use_rep=emb_key, key_added=f"neighbors_{method_name}")
    sc.tl.umap(adata, neighbors_key=f"neighbors_{method_name}")

    color_keys_present = [k for k in color_keys if k in adata.obs]
    n_colors = len(color_keys_present)
    if n_colors == 0:
        return

    fig, axes = plt.subplots(1, n_colors, figsize=(6 * n_colors, 5))
    if n_colors == 1:
        axes = [axes]
    fig.suptitle(f"{milestone}: {method_name}", fontsize=14, fontweight="bold")
    for j, key in enumerate(color_keys_present):
        sc.pl.umap(adata, color=key, ax=axes[j], show=False, title=key)
    plt.tight_layout()
    fname = f"{milestone}_umap_{method_name.lower().replace(' ', '_').replace('+', '_')}.png"
    fig.savefig(figdir / fname, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {fname}")


# ============================================================
# Process M1 and M2 (M0 has no cell type labels)
# ============================================================

results_all = []
color_keys = ["platform", "panel_variant", "patient_id", "disease", "celltype_broad"]

for milestone, config in [
    ("m1", {
        "adata_path": WORKDIR / "results" / "m1_benchmark" / "adata.m1.h5ad",
        "batch_key": "platform",
        "label_key": "celltype_broad",
        "n_latent": 10,
        "max_epochs": 100,
        "figdir": WORKDIR / "figures" / "m1",
    }),
    ("m2", {
        "adata_path": WORKDIR / "results" / "m2_benchmark" / "adata.m2.h5ad",
        "batch_key": "platform",
        "label_key": "celltype_broad",
        "n_latent": 20,
        "max_epochs": 100,
        "figdir": WORKDIR / "figures" / "m2",
    }),
]:
    print(f"\n{'='*60}")
    print(f"  {milestone.upper()}: resolVI semi-supervised")
    print(f"{'='*60}")

    adata_path = config["adata_path"]
    if not adata_path.exists():
        print(f"  SKIP: {adata_path} not found")
        continue

    adata = ad.read_h5ad(adata_path)
    print(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"  Embeddings: {list(adata.obsm.keys())}")

    batch_key = config["batch_key"]
    label_key = config["label_key"]
    figdir = config["figdir"]
    figdir.mkdir(parents=True, exist_ok=True)

    # Skip if already computed
    if "X_resolvi_ss" in adata.obsm:
        print(f"  resolVI semi-supervised already exists, computing metrics only")
        batch_asw, ct_asw, ct_broad_asw = compute_metrics(
            adata, "X_resolvi_ss", batch_key, "celltype"
        )
        batch_score = 1 - abs(batch_asw)
        entropy, n_clusters = compute_platform_entropy(adata, "X_resolvi_ss", batch_key)
    else:
        try:
            model, adata_rv = run_resolvi_semisup(
                adata, batch_key, label_key,
                n_latent=config["n_latent"],
                max_epochs=config["max_epochs"],
            )

            # Get latent representation
            latent = model.get_latent_representation()
            print(f"  resolVI-SS latent: {latent.shape}")
            adata.obsm["X_resolvi_ss"] = latent.astype(np.float32)

            # Get cell type predictions
            try:
                predictions = model.predict(adata_rv)
                if isinstance(predictions, pd.DataFrame):
                    adata.obs["resolvi_ss_pred"] = predictions.idxmax(axis=1).values
                    adata.obsm["resolvi_ss_probs"] = predictions.values

                    # Prediction accuracy
                    orig = adata.obs[label_key].astype(str)
                    pred = adata.obs["resolvi_ss_pred"].astype(str)
                    valid = (orig != "nan") & (orig != "Unknown")
                    if valid.sum() > 0:
                        acc = (orig[valid] == pred[valid]).mean()
                        print(f"  Prediction accuracy: {acc:.3f} ({valid.sum()} cells)")
                        for platform in adata.obs[batch_key].unique():
                            pmask = valid & (adata.obs[batch_key] == platform)
                            if pmask.sum() > 0:
                                pacc = (orig[pmask] == pred[pmask]).mean()
                                print(f"    {platform}: {pacc:.3f} ({pmask.sum()} cells)")
                else:
                    print(f"  predict() returned {type(predictions)}, storing as-is")
                    adata.obs["resolvi_ss_pred"] = predictions
            except Exception as e:
                print(f"  predict() failed (non-fatal): {e}")

            # Metrics
            batch_asw, ct_asw, ct_broad_asw = compute_metrics(
                adata, "X_resolvi_ss", batch_key, "celltype"
            )
            batch_score = 1 - abs(batch_asw)
            entropy, n_clusters = compute_platform_entropy(adata, "X_resolvi_ss", batch_key)

        except Exception as e:
            print(f"  resolVI semi-supervised FAILED: {e}")
            import traceback
            traceback.print_exc()
            continue

    print(f"  resolVI-SS: batch_ASW={batch_asw:.3f}, batch_score={batch_score:.3f}, "
          f"ct_broad_ASW={ct_broad_asw if ct_broad_asw is not None else 'N/A'}, "
          f"entropy={entropy:.3f}")

    results_all.append({
        "milestone": milestone, "method": "resolVI-SS",
        "batch_asw": batch_asw, "celltype_asw": ct_asw,
        "celltype_broad_asw": ct_broad_asw, "batch_score": batch_score,
        "platform_entropy": entropy,
    })

    # UMAP
    plot_umaps(adata, "X_resolvi_ss", "resolVI_SS", milestone, figdir, batch_key, color_keys)

    # --- Compare against existing methods ---
    print(f"\n  --- Comparison for {milestone.upper()} ---")
    comparison = []
    for method, key in [
        ("resolVI (unsup)", "X_resolvi"),
        ("resolVI-SS", "X_resolvi_ss"),
        ("scANVI", "X_scanvi"),
        ("scVI", "X_scvi"),
        ("Harmony", "X_harmony"),
    ]:
        if key in adata.obsm:
            b_asw, c_asw, cb_asw = compute_metrics(adata, key, batch_key, "celltype")
            comparison.append({
                "method": method, "batch_score": 1 - abs(b_asw),
                "ct_broad_asw": cb_asw,
            })
    if comparison:
        comp_df = pd.DataFrame(comparison).sort_values("ct_broad_asw", ascending=False)
        print(comp_df.to_string(index=False))

    # Save adata with new embedding
    print(f"\n  Saving {milestone} adata...")
    adata.write_h5ad(adata_path)
    print(f"  Saved: {adata_path} ({adata.n_obs} cells)")

# ============================================================
# Summary
# ============================================================
print(f"\n{'='*60}")
print("RESOLVI SEMI-SUPERVISED COMPLETE")
print(f"{'='*60}")

if results_all:
    results_df = pd.DataFrame(results_all)
    print(results_df.to_string(index=False))

    out_path = WORKDIR / "results" / "resolvi_semisup_benchmark.csv"
    results_df.to_csv(out_path, index=False)
    print(f"\nSaved: {out_path}")
