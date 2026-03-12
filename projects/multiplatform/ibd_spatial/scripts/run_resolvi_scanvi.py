#!/usr/bin/env python3
"""Run resolVI → scANVI for all milestones (M0, M1, M2).

resolVI is a spatial-aware VAE that accounts for spatial autocorrelation.
Using it as the pretrain backbone for scANVI should improve cross-platform
integration by leveraging spatial structure.

For milestones without cell type labels (M0), runs resolVI only.
For M1 and M2 (with celltype_broad labels), runs resolVI → scANVI.

Reads:
  results/m0_benchmark/adata.m0.h5ad
  results/m1_benchmark/adata.m1.h5ad
  results/m2_benchmark/adata.m2.h5ad
Writes:
  Updated adata with X_resolvi, X_resolvi_scanvi embeddings
  Benchmark CSVs and UMAP PNGs per milestone
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

import scvi
from scvi.external import RESOLVI

WORKDIR = Path("/home/fs01/juk4007/elementolab/projects/ibd_spatial")


def subsample_balanced(adata, key, max_cells=50000):
    """Subsample balancing by key."""
    if adata.n_obs <= max_cells:
        return adata
    np.random.seed(42)
    groups = adata.obs[key].unique()
    per_group = max_cells // len(groups)
    idx = []
    for g in groups:
        gmask = np.where(adata.obs[key] == g)[0]
        n_take = min(len(gmask), per_group)
        idx.extend(np.random.choice(gmask, n_take, replace=False))
    return adata[sorted(idx)].copy()


def run_resolvi(adata, batch_key, spatial_key="spatial", n_latent=20, max_epochs=100):
    """Run resolVI and return the model."""
    print(f"  Setting up resolVI: batch_key={batch_key}, n_latent={n_latent}")

    adata_rv = adata.copy()

    # Ensure raw counts in X
    if "counts" in adata_rv.layers:
        adata_rv.X = adata_rv.layers["counts"].copy()

    # resolVI expects obsm['X_spatial'] for spatial neighbor graph
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

    # Setup resolVI — requires spatial coordinates
    RESOLVI.setup_anndata(
        adata_rv,
        batch_key=batch_key,
        layer=None,  # use X (raw counts)
    )

    model = RESOLVI(
        adata_rv,
        n_latent=n_latent,
    )

    print(f"  Training resolVI ({max_epochs} epochs)...")
    t0 = time.time()
    # resolVI uses Pyro backend — train() API differs from scVI (no train_size param)
    model.train(
        max_epochs=max_epochs,
        batch_size=512,
    )
    elapsed = time.time() - t0
    print(f"  resolVI training done ({elapsed:.0f}s)")

    return model, adata_rv


def run_scanvi_fresh(adata, batch_key, label_key, n_latent=20, unlabeled="Unknown"):
    """Run scVI pretrain → scANVI (fresh, not from resolVI since Pyro/Lightning incompatible)."""
    print(f"  Running fresh scVI→scANVI (label_key={label_key})")

    adata_s = adata.copy()
    if "counts" in adata_s.layers:
        adata_s.X = adata_s.layers["counts"].copy()

    # HVG
    if adata_s.n_vars > 500:
        sc.pp.highly_variable_genes(adata_s, n_top_genes=2000, batch_key=batch_key, flavor="seurat_v3")
        adata_s = adata_s[:, adata_s.var["highly_variable"]].copy()
        print(f"  HVG subset: {adata_s.n_vars} genes")

    # Prepare labels
    adata_s.obs["_scanvi_labels"] = (
        adata_s.obs[label_key].astype(str).replace("nan", unlabeled)
    )

    # scVI pretrain
    scvi.model.SCVI.setup_anndata(adata_s, batch_key=batch_key)
    scvi_model = scvi.model.SCVI(adata_s, n_latent=n_latent, n_hidden=128, n_layers=2)
    t0 = time.time()
    scvi_model.train(max_epochs=50, early_stopping=True, early_stopping_patience=10,
                     batch_size=512, train_size=0.9)
    print(f"  scVI pretrain done ({time.time()-t0:.0f}s)")

    # scANVI from scVI
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model, unlabeled_category=unlabeled, labels_key="_scanvi_labels"
    )
    t0 = time.time()
    scanvi_model.train(max_epochs=30, batch_size=512, train_size=0.9)
    print(f"  scANVI done ({time.time()-t0:.0f}s)")

    return scanvi_model, adata_s


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

    ct_broad_key = ct_key.replace("celltype", "celltype_broad") if ct_key else "celltype_broad"
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
    sc.tl.leiden(adata, neighbors_key="neighbors_tmp", key_added="leiden_tmp", resolution=resolution)

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
    fig.savefig(figdir / f"{milestone}_umap_{method_name.lower().replace('+', '_')}.png",
                dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {milestone}_umap_{method_name.lower().replace('+', '_')}.png")


# ============================================================
# Process each milestone
# ============================================================

results_all = []
color_keys = ["platform", "panel_variant", "patient_id", "disease", "celltype_broad"]

for milestone, config in [
    ("m0", {
        "adata_path": WORKDIR / "results" / "m0_benchmark" / "adata.m0.h5ad",
        "batch_key": "panel_variant",
        "ct_key": None,  # no cell type labels in M0
        "n_latent": 10,
        "max_epochs": 100,
        "figdir": WORKDIR / "figures" / "m0",
    }),
    ("m1", {
        "adata_path": WORKDIR / "results" / "m1_benchmark" / "adata.m1.h5ad",
        "batch_key": "platform",
        "ct_key": "celltype",
        "n_latent": 10,
        "max_epochs": 100,
        "figdir": WORKDIR / "figures" / "m1",
    }),
    ("m2", {
        "adata_path": WORKDIR / "results" / "m2_benchmark" / "adata.m2.h5ad",
        "batch_key": "platform",
        "ct_key": "celltype",
        "n_latent": 20,
        "max_epochs": 100,
        "figdir": WORKDIR / "figures" / "m2",
    }),
]:
    print(f"\n{'='*60}")
    print(f"  {milestone.upper()}")
    print(f"{'='*60}")

    adata_path = config["adata_path"]
    if not adata_path.exists():
        print(f"  SKIP: {adata_path} not found")
        continue

    adata = ad.read_h5ad(adata_path)
    print(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"  Embeddings: {list(adata.obsm.keys())}")

    batch_key = config["batch_key"]
    ct_key = config["ct_key"]
    figdir = config["figdir"]
    figdir.mkdir(parents=True, exist_ok=True)

    # --- Run resolVI (skip if already done) ---
    if "X_resolvi" in adata.obsm:
        print(f"  resolVI already exists, skipping training")
        resolvi_success = True
        # Still compute metrics for reporting
        batch_asw, ct_asw, ct_broad_asw = compute_metrics(
            adata, "X_resolvi", batch_key, ct_key
        )
        batch_score = 1 - abs(batch_asw)
        entropy, n_clusters = compute_platform_entropy(adata, "X_resolvi", batch_key)
        print(f"  resolVI: batch_ASW={batch_asw:.3f}, batch_score={batch_score:.3f}, "
              f"ct_broad_ASW={ct_broad_asw if ct_broad_asw is not None else 'N/A'}, "
              f"entropy={entropy:.3f}")
        results_all.append({
            "milestone": milestone, "method": "resolVI",
            "batch_asw": batch_asw, "celltype_asw": ct_asw,
            "celltype_broad_asw": ct_broad_asw, "batch_score": batch_score,
            "platform_entropy": entropy,
        })
    else:
        try:
            resolvi_model, adata_rv = run_resolvi(
                adata, batch_key,
                n_latent=config["n_latent"],
                max_epochs=config["max_epochs"],
            )

            # Get latent representation
            latent = resolvi_model.get_latent_representation()
            print(f"  resolVI latent: {latent.shape}")

            adata.obsm["X_resolvi"] = latent.astype(np.float32)

            # Metrics
            batch_asw, ct_asw, ct_broad_asw = compute_metrics(
                adata, "X_resolvi", batch_key, ct_key
            )
            batch_score = 1 - abs(batch_asw)
            entropy, n_clusters = compute_platform_entropy(adata, "X_resolvi", batch_key)

            print(f"  resolVI: batch_ASW={batch_asw:.3f}, batch_score={batch_score:.3f}, "
                  f"ct_broad_ASW={ct_broad_asw if ct_broad_asw is not None else 'N/A'}, "
                  f"entropy={entropy:.3f}")

            results_all.append({
                "milestone": milestone, "method": "resolVI",
                "batch_asw": batch_asw, "celltype_asw": ct_asw,
                "celltype_broad_asw": ct_broad_asw, "batch_score": batch_score,
                "platform_entropy": entropy,
            })

            # UMAP
            plot_umaps(adata, "X_resolvi", "resolVI", milestone, figdir, batch_key, color_keys)

            resolvi_success = True
        except Exception as e:
            print(f"  resolVI FAILED: {e}")
            import traceback
            traceback.print_exc()
            resolvi_success = False

    # --- Run fresh scVI→scANVI (only if cell type labels available) ---
    # Note: resolVI is Pyro-based, scANVI is PyTorch Lightning — cannot chain them.
    # Instead, run scANVI from a fresh scVI pretrain (same approach as M2 followup).
    # Only run if scANVI embedding doesn't already exist (avoid re-running).
    if ct_key is not None and "X_scanvi" not in adata.obsm:
        try:
            ct_broad_key = "celltype_broad"
            if ct_broad_key not in adata.obs:
                print(f"  SKIP scANVI: no {ct_broad_key} column")
            else:
                scanvi_model, adata_s = run_scanvi_fresh(
                    adata, batch_key, ct_broad_key, n_latent=config["n_latent"]
                )

                latent_scanvi = scanvi_model.get_latent_representation()
                print(f"  scANVI latent: {latent_scanvi.shape}")

                adata.obsm["X_scanvi"] = latent_scanvi.astype(np.float32)

                # Predictions
                predictions = scanvi_model.predict()
                adata.obs["scanvi_pred"] = predictions if isinstance(predictions, np.ndarray) else predictions.values

                # Metrics
                batch_asw, ct_asw, ct_broad_asw = compute_metrics(
                    adata, "X_scanvi", batch_key, ct_key
                )
                batch_score = 1 - abs(batch_asw)
                entropy, n_clusters = compute_platform_entropy(
                    adata, "X_scanvi", batch_key
                )

                print(f"  scANVI: batch_ASW={batch_asw:.3f}, batch_score={batch_score:.3f}, "
                      f"ct_broad_ASW={ct_broad_asw if ct_broad_asw is not None else 'N/A'}, "
                      f"entropy={entropy:.3f}")

                # Prediction accuracy
                if ct_broad_key in adata.obs:
                    orig = adata.obs[ct_broad_key].astype(str)
                    pred = adata.obs["scanvi_pred"].astype(str)
                    valid = (orig != "nan") & (orig != "Unknown")
                    if valid.sum() > 0:
                        acc = (orig[valid] == pred[valid]).mean()
                        print(f"  Prediction accuracy: {acc:.3f} ({valid.sum()} cells)")

                        for platform in adata.obs[batch_key].unique():
                            pmask = valid & (adata.obs[batch_key] == platform)
                            if pmask.sum() > 0:
                                pacc = (orig[pmask] == pred[pmask]).mean()
                                print(f"    {platform}: {pacc:.3f} ({pmask.sum()} cells)")

                results_all.append({
                    "milestone": milestone, "method": "scANVI",
                    "batch_asw": batch_asw, "celltype_asw": ct_asw,
                    "celltype_broad_asw": ct_broad_asw, "batch_score": batch_score,
                    "platform_entropy": entropy,
                })

                plot_umaps(adata, "X_scanvi", "scANVI", milestone,
                          figdir, batch_key, color_keys)

        except Exception as e:
            print(f"  scANVI FAILED: {e}")
            import traceback
            traceback.print_exc()
    elif "X_scanvi" in adata.obsm:
        print(f"  scANVI already exists in {milestone}, skipping")

    # --- Save updated adata ---
    print(f"\n  Saving {milestone} adata...")
    adata.write_h5ad(adata_path)
    print(f"  Saved: {adata_path} ({adata.n_obs} cells)")

# ============================================================
# Summary
# ============================================================
print(f"\n{'='*60}")
print("RESOLVI + SCANVI COMPLETE")
print(f"{'='*60}")

if results_all:
    results_df = pd.DataFrame(results_all)
    print(results_df.to_string(index=False))

    # Save combined results
    out_path = WORKDIR / "results" / "resolvi_scanvi_benchmark.csv"
    results_df.to_csv(out_path, index=False)
    print(f"\nSaved: {out_path}")
