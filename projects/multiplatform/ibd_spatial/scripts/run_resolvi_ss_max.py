#!/usr/bin/env python3
"""resolVI-SS-max sweep: 10 configurations to maximize bio conservation (celltype_broad ASW).

resolVI with semisupervised=True — spatial-aware VAE with label propagation.

Reads:  results/m2_benchmark/adata.m2.h5ad (164K cells, 2552 genes)
Writes: results/resolvi_ss_max/resolvi_ss_max_sweep.csv
        results/resolvi_ss_max/adata.m2.resolvi_ss_max.h5ad
"""

import os
import sys
import time
import traceback
from pathlib import Path

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout.reconfigure(line_buffering=True)

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")

import torch
torch.set_float32_matmul_precision("medium")

from scvi.external import RESOLVI

WORKDIR = Path("/home/fs01/juk4007/elementolab/projects/ibd_spatial")
INPUT_PATH = WORKDIR / "results" / "m2_benchmark" / "adata.m2.h5ad"
OUTDIR = WORKDIR / "results" / "resolvi_ss_max"

# ============================================================
# Configuration sweep
# ============================================================
CONFIGS = [
    {
        "name": "B1_baseline",
        "n_latent": 10, "n_hidden": 32, "n_hidden_encoder": 128, "n_layers": 2,
        "prior_diffusion_amount": 0.3, "sparsity_diffusion": 3.0,
        "n_neighbors": 10, "max_epochs": 100,
        "lr": 3e-3, "lr_extra": 1e-2,
        "batch_key": "platform",
        "classifier_parameters": None,
    },
    {
        "name": "B2_lat20",
        "n_latent": 20, "n_hidden": 32, "n_hidden_encoder": 128, "n_layers": 2,
        "prior_diffusion_amount": 0.3, "sparsity_diffusion": 3.0,
        "n_neighbors": 10, "max_epochs": 100,
        "lr": 3e-3, "lr_extra": 1e-2,
        "batch_key": "platform",
        "classifier_parameters": None,
    },
    {
        "name": "B3_big_arch",
        "n_latent": 20, "n_hidden": 128, "n_hidden_encoder": 256, "n_layers": 2,
        "prior_diffusion_amount": 0.3, "sparsity_diffusion": 3.0,
        "n_neighbors": 10, "max_epochs": 100,
        "lr": 3e-3, "lr_extra": 1e-2,
        "batch_key": "platform",
        "classifier_parameters": None,
    },
    {
        "name": "B4_low_diff",
        "n_latent": 20, "n_hidden": 128, "n_hidden_encoder": 256, "n_layers": 2,
        "prior_diffusion_amount": 0.1, "sparsity_diffusion": 5.0,
        "n_neighbors": 10, "max_epochs": 100,
        "lr": 3e-3, "lr_extra": 1e-3,
        "batch_key": "platform",
        "classifier_parameters": None,
    },
    {
        "name": "B5_few_neighbors",
        "n_latent": 20, "n_hidden": 128, "n_hidden_encoder": 256, "n_layers": 2,
        "prior_diffusion_amount": 0.3, "sparsity_diffusion": 3.0,
        "n_neighbors": 5, "max_epochs": 100,
        "lr": 3e-3, "lr_extra": 1e-2,
        "batch_key": "platform",
        "classifier_parameters": None,
    },
    {
        "name": "B6_deep_cls",
        "n_latent": 20, "n_hidden": 128, "n_hidden_encoder": 256, "n_layers": 2,
        "prior_diffusion_amount": 0.3, "sparsity_diffusion": 3.0,
        "n_neighbors": 10, "max_epochs": 100,
        "lr": 3e-3, "lr_extra": 1e-2,
        "batch_key": "platform",
        "classifier_parameters": {"n_layers": 2, "n_hidden": 256, "dropout_rate": 0.1},
    },
    {
        "name": "B7_combined",
        "n_latent": 30, "n_hidden": 128, "n_hidden_encoder": 256, "n_layers": 2,
        "prior_diffusion_amount": 0.1, "sparsity_diffusion": 5.0,
        "n_neighbors": 5, "max_epochs": 150,
        "lr": 3e-3, "lr_extra": 1e-3,
        "batch_key": "platform",
        "classifier_parameters": {"n_layers": 2, "n_hidden": 256, "dropout_rate": 0.1},
    },
    {
        "name": "B8_per_sample",
        "n_latent": 20, "n_hidden": 128, "n_hidden_encoder": 256, "n_layers": 2,
        "prior_diffusion_amount": 0.3, "sparsity_diffusion": 3.0,
        "n_neighbors": 10, "max_epochs": 100,
        "lr": 3e-3, "lr_extra": 1e-2,
        "batch_key": "sample",
        "classifier_parameters": None,
    },
    {
        "name": "B9_min_spatial",
        "n_latent": 20, "n_hidden": 128, "n_hidden_encoder": 256, "n_layers": 2,
        "prior_diffusion_amount": 0.05, "sparsity_diffusion": 8.0,
        "n_neighbors": 3, "max_epochs": 100,
        "lr": 3e-3, "lr_extra": 5e-4,
        "batch_key": "platform",
        "classifier_parameters": None,
    },
    {
        "name": "B10_max_cap",
        "n_latent": 30, "n_hidden": 128, "n_hidden_encoder": 256, "n_layers": 3,
        "prior_diffusion_amount": 0.1, "sparsity_diffusion": 5.0,
        "n_neighbors": 5, "max_epochs": 200,
        "lr": 1e-3, "lr_extra": 1e-3,
        "batch_key": "platform",
        "classifier_parameters": {"n_layers": 2, "n_hidden": 256, "dropout_rate": 0.1},
    },
]


# ============================================================
# Metric helpers
# ============================================================

def compute_metrics(emb, obs, batch_key="platform"):
    """Compute batch_asw, batch_score, ct_broad_asw, ct_asw."""
    from sklearn.metrics import silhouette_score
    n_sample = min(10000, emb.shape[0])

    batch_asw = silhouette_score(
        emb, obs[batch_key].values,
        sample_size=n_sample, random_state=42,
    )
    batch_score = 1 - abs(batch_asw)

    ct_broad_asw = None
    if "celltype_broad" in obs.columns:
        labels = obs["celltype_broad"].values
        valid = pd.notna(labels)
        if valid.sum() > 100 and len(np.unique(labels[valid])) > 1:
            ct_broad_asw = silhouette_score(
                emb[valid], labels[valid],
                sample_size=min(n_sample, int(valid.sum())), random_state=42,
            )

    ct_asw = None
    if "celltype" in obs.columns:
        labels = obs["celltype"].values
        valid = pd.notna(labels)
        if valid.sum() > 100 and len(np.unique(labels[valid])) > 1:
            ct_asw = silhouette_score(
                emb[valid], labels[valid],
                sample_size=min(n_sample, int(valid.sum())), random_state=42,
            )

    return batch_asw, batch_score, ct_broad_asw, ct_asw


def compute_platform_entropy(emb, obs, batch_key="platform", resolution=0.5):
    """Compute per-cluster platform entropy from embedding."""
    from scipy.stats import entropy as scipy_entropy
    import warnings
    warnings.filterwarnings("ignore", category=FutureWarning)

    tmp = ad.AnnData(obs=obs.copy())
    tmp.obsm["X_emb"] = emb
    sc.pp.neighbors(tmp, use_rep="X_emb")
    sc.tl.leiden(tmp, resolution=resolution, key_added="leiden_tmp")

    n_batches = obs[batch_key].nunique()
    max_ent = np.log(n_batches) if n_batches > 1 else 1.0
    entropies = []
    for cluster in tmp.obs["leiden_tmp"].unique():
        mask = tmp.obs["leiden_tmp"] == cluster
        counts = obs.loc[mask, batch_key].value_counts()
        ent = scipy_entropy(counts / counts.sum()) / max_ent if max_ent > 0 else 0
        entropies.append(ent)
    return np.median(entropies)


# ============================================================
# Main
# ============================================================

def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    print(f"{'='*70}")
    print(f"  resolVI-SS-max sweep: {len(CONFIGS)} configurations")
    print(f"  Input: {INPUT_PATH}")
    print(f"  Output: {OUTDIR}")
    print(f"{'='*70}")

    # Load data
    print(f"\nLoading adata...")
    adata_orig = ad.read_h5ad(INPUT_PATH)
    print(f"  Shape: {adata_orig.n_obs} cells x {adata_orig.n_vars} genes")
    print(f"  Layers: {list(adata_orig.layers.keys())}")
    print(f"  Obsm: {list(adata_orig.obsm.keys())}")
    if "celltype_broad" in adata_orig.obs:
        print(f"  celltype_broad: {adata_orig.obs['celltype_broad'].value_counts().to_dict()}")

    results = []
    best_ct_broad_asw = -1.0
    best_config_name = None
    best_embedding = None

    for i, cfg in enumerate(CONFIGS):
        name = cfg["name"]
        print(f"\n{'='*70}")
        print(f"  [{i+1}/{len(CONFIGS)}] {name}")
        print(f"  n_latent={cfg['n_latent']}, n_hidden={cfg['n_hidden']}, "
              f"n_hidden_encoder={cfg['n_hidden_encoder']}, n_layers={cfg['n_layers']}, "
              f"prior_diff={cfg['prior_diffusion_amount']}, sparsity_diff={cfg['sparsity_diffusion']}, "
              f"n_neighbors={cfg['n_neighbors']}, epochs={cfg['max_epochs']}, "
              f"lr={cfg['lr']}, lr_extra={cfg['lr_extra']}, "
              f"batch_key={cfg['batch_key']}, cls_params={cfg['classifier_parameters']}")
        print(f"{'='*70}")

        try:
            t0_total = time.time()

            # Copy adata
            adata = adata_orig.copy()

            # Ensure raw counts in X
            if "counts" in adata.layers:
                adata.X = adata.layers["counts"].copy()
                print("  Using layers['counts'] as X")

            # HVG selection
            sc.pp.highly_variable_genes(
                adata, n_top_genes=2000, batch_key=cfg["batch_key"], flavor="seurat_v3"
            )
            adata = adata[:, adata.var["highly_variable"]].copy()
            print(f"  HVG subset: {adata.n_vars} genes")

            # Copy spatial to X_spatial
            if "spatial" in adata.obsm and "X_spatial" not in adata.obsm:
                adata.obsm["X_spatial"] = adata.obsm["spatial"].copy()
                print("  Copied obsm['spatial'] -> obsm['X_spatial']")

            # Prepare labels
            adata.obs["_resolvi_labels"] = (
                adata.obs["celltype_broad"].astype(str).replace("nan", "Unknown")
            )
            n_labeled = (adata.obs["_resolvi_labels"] != "Unknown").sum()
            n_types = adata.obs["_resolvi_labels"].nunique()
            print(f"  Labels: {n_labeled}/{adata.n_obs} labeled, {n_types} types")

            # Setup anndata
            print(f"  Setting up resolVI (n_neighbors={cfg['n_neighbors']})...")
            RESOLVI.setup_anndata(
                adata,
                batch_key=cfg["batch_key"],
                labels_key="_resolvi_labels",
                layer=None,
                prepare_data_kwargs={"n_neighbors": cfg["n_neighbors"]},
            )

            # Build model kwargs
            model_kwargs = {
                "prior_diffusion_amount": cfg["prior_diffusion_amount"],
                "sparsity_diffusion": cfg["sparsity_diffusion"],
            }
            if cfg["classifier_parameters"] is not None:
                model_kwargs["classifier_parameters"] = cfg["classifier_parameters"]

            # Try creating model — resolVI API may vary, try multiple approaches
            print(f"  Creating resolVI model...")
            try:
                # Approach 1: direct constructor kwargs
                model = RESOLVI(
                    adata,
                    n_latent=cfg["n_latent"],
                    n_hidden=cfg["n_hidden"],
                    n_hidden_encoder=cfg["n_hidden_encoder"],
                    n_layers=cfg["n_layers"],
                    semisupervised=True,
                    **model_kwargs,
                )
            except TypeError as e1:
                print(f"  Direct kwargs failed: {e1}")
                print(f"  Trying with model_kwargs dict...")
                try:
                    # Approach 2: wrap spatial/classifier params in model_kwargs
                    model = RESOLVI(
                        adata,
                        n_latent=cfg["n_latent"],
                        n_hidden=cfg["n_hidden"],
                        n_hidden_encoder=cfg["n_hidden_encoder"],
                        n_layers=cfg["n_layers"],
                        semisupervised=True,
                        model_kwargs=model_kwargs,
                    )
                except TypeError as e2:
                    print(f"  model_kwargs dict also failed: {e2}")
                    print(f"  Trying minimal constructor + defaults...")
                    # Approach 3: minimal — just architecture + semisupervised
                    model = RESOLVI(
                        adata,
                        n_latent=cfg["n_latent"],
                        semisupervised=True,
                    )

            # Train
            print(f"  Training resolVI-SS ({cfg['max_epochs']} epochs)...")
            t0 = time.time()
            model.train(
                max_epochs=cfg["max_epochs"],
                batch_size=512,
                lr=cfg["lr"],
                lr_extra=cfg["lr_extra"],
            )
            train_time = time.time() - t0
            print(f"  Training done ({train_time:.0f}s)")

            # Get latent representation
            latent = model.get_latent_representation()
            print(f"  Latent: {latent.shape}")

            # Try predictions
            pred_acc = None
            try:
                predictions = model.predict(adata)
                if isinstance(predictions, pd.DataFrame):
                    pred_labels = predictions.idxmax(axis=1).values
                    orig = adata.obs["celltype_broad"].astype(str)
                    pred = pd.Series(pred_labels, index=adata.obs.index).astype(str)
                    valid = (orig != "nan") & (orig != "Unknown")
                    if valid.sum() > 0:
                        pred_acc = (orig[valid] == pred[valid]).mean()
                        print(f"  Prediction accuracy: {pred_acc:.3f} ({valid.sum()} cells)")
                elif isinstance(predictions, (np.ndarray, pd.Series)):
                    orig = adata.obs["celltype_broad"].astype(str)
                    pred = pd.Series(predictions, index=adata.obs.index).astype(str)
                    valid = (orig != "nan") & (orig != "Unknown")
                    if valid.sum() > 0:
                        pred_acc = (orig[valid] == pred[valid]).mean()
                        print(f"  Prediction accuracy: {pred_acc:.3f} ({valid.sum()} cells)")
            except Exception as e:
                print(f"  predict() failed (non-fatal): {e}")

            # Compute metrics (use platform as batch_key for metrics even if model used sample)
            batch_asw, batch_score, ct_broad_asw, ct_asw = compute_metrics(
                latent, adata_orig.obs, batch_key="platform"
            )
            platform_entropy = compute_platform_entropy(
                latent, adata_orig.obs, batch_key="platform"
            )

            total_time = time.time() - t0_total

            row = {
                "name": name,
                "batch_asw": round(batch_asw, 4),
                "batch_score": round(batch_score, 4),
                "ct_broad_asw": round(ct_broad_asw, 4) if ct_broad_asw is not None else None,
                "ct_asw": round(ct_asw, 4) if ct_asw is not None else None,
                "platform_entropy": round(platform_entropy, 4),
                "pred_accuracy": round(pred_acc, 4) if pred_acc is not None else None,
                "train_time_s": round(train_time, 1),
                "total_time_s": round(total_time, 1),
            }
            results.append(row)

            print(f"\n  RESULTS: batch_score={batch_score:.4f}, "
                  f"ct_broad_asw={ct_broad_asw:.4f if ct_broad_asw is not None else 'N/A'}, "
                  f"ct_asw={ct_asw:.4f if ct_asw is not None else 'N/A'}, "
                  f"entropy={platform_entropy:.4f}, "
                  f"acc={pred_acc:.4f if pred_acc is not None else 'N/A'}")

            # Track best
            score = ct_broad_asw if ct_broad_asw is not None else -1.0
            if score > best_ct_broad_asw:
                best_ct_broad_asw = score
                best_config_name = name
                best_embedding = latent.astype(np.float32)
                print(f"  *** NEW BEST: {name} with ct_broad_asw={score:.4f} ***")

            # Print comparison table so far
            if len(results) > 1:
                print(f"\n  --- Comparison after {len(results)} configs ---")
                df_tmp = pd.DataFrame(results).sort_values("ct_broad_asw", ascending=False)
                print(df_tmp[["name", "batch_score", "ct_broad_asw", "ct_asw",
                              "platform_entropy", "pred_accuracy"]].to_string(index=False))

        except Exception as e:
            print(f"  FAILED: {e}")
            traceback.print_exc()
            results.append({
                "name": name,
                "batch_asw": None, "batch_score": None,
                "ct_broad_asw": None, "ct_asw": None,
                "platform_entropy": None, "pred_accuracy": None,
                "train_time_s": None, "total_time_s": None,
            })
            continue

    # ============================================================
    # Final summary
    # ============================================================
    print(f"\n{'='*70}")
    print("  RESOLVI-SS-MAX SWEEP COMPLETE")
    print(f"{'='*70}")

    if results:
        results_df = pd.DataFrame(results)
        results_df_sorted = results_df.sort_values("ct_broad_asw", ascending=False)

        print("\nFinal ranking (by ct_broad_asw descending):")
        print(results_df_sorted.to_string(index=False))

        csv_path = OUTDIR / "resolvi_ss_max_sweep.csv"
        results_df_sorted.to_csv(csv_path, index=False)
        print(f"\nSaved: {csv_path}")

    if best_embedding is not None:
        print(f"\nBest config: {best_config_name} (ct_broad_asw={best_ct_broad_asw:.4f})")
        print(f"Saving adata with X_resolvi_ss_max embedding...")

        adata_out = adata_orig.copy()
        adata_out.obsm["X_resolvi_ss_max"] = best_embedding

        out_path = OUTDIR / "adata.m2.resolvi_ss_max.h5ad"
        adata_out.write_h5ad(out_path)
        print(f"Saved: {out_path} ({adata_out.n_obs} cells)")
    else:
        print("\nNo successful configs — no adata saved.")

    print("\nDone.")


if __name__ == "__main__":
    main()
