#!/usr/bin/env python3
"""scANVI-max sweep: 12 configurations to maximize bio conservation (celltype_broad ASW).

Reads:  results/m2_benchmark/adata.m2.h5ad (164K cells, 2552 genes)
Writes: results/scanvi_max/scanvi_max_sweep.csv
        results/scanvi_max/adata.m2.scanvi_max.h5ad
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

import scvi

WORKDIR = Path("/home/fs01/juk4007/elementolab/projects/ibd_spatial")
INPUT_PATH = WORKDIR / "results" / "m2_benchmark" / "adata.m2.h5ad"
OUTDIR = WORKDIR / "results" / "scanvi_max"

# ============================================================
# Configuration sweep
# ============================================================
CONFIGS = [
    {
        "name": "A1_baseline",
        "n_latent": 20, "n_hidden": 128, "n_layers": 2,
        "gene_likelihood": "zinb", "dispersion": "gene",
        "classification_ratio": 50,
        "pretrain_epochs": 50, "finetune_epochs": 30,
        "label_key": "celltype_broad",
        "n_samples_per_label": None,
        "plan_kwargs_extra": {},
    },
    {
        "name": "A2_cls100",
        "n_latent": 20, "n_hidden": 128, "n_layers": 2,
        "gene_likelihood": "zinb", "dispersion": "gene",
        "classification_ratio": 100,
        "pretrain_epochs": 50, "finetune_epochs": 30,
        "label_key": "celltype_broad",
        "n_samples_per_label": None,
        "plan_kwargs_extra": {},
    },
    {
        "name": "A3_cls200",
        "n_latent": 20, "n_hidden": 128, "n_layers": 2,
        "gene_likelihood": "zinb", "dispersion": "gene",
        "classification_ratio": 200,
        "pretrain_epochs": 50, "finetune_epochs": 30,
        "label_key": "celltype_broad",
        "n_samples_per_label": None,
        "plan_kwargs_extra": {},
    },
    {
        "name": "A4_big_arch",
        "n_latent": 30, "n_hidden": 256, "n_layers": 2,
        "gene_likelihood": "zinb", "dispersion": "gene",
        "classification_ratio": 100,
        "pretrain_epochs": 100, "finetune_epochs": 50,
        "label_key": "celltype_broad",
        "n_samples_per_label": None,
        "plan_kwargs_extra": {},
    },
    {
        "name": "A5_nb",
        "n_latent": 30, "n_hidden": 256, "n_layers": 2,
        "gene_likelihood": "nb", "dispersion": "gene",
        "classification_ratio": 100,
        "pretrain_epochs": 100, "finetune_epochs": 50,
        "label_key": "celltype_broad",
        "n_samples_per_label": None,
        "plan_kwargs_extra": {},
    },
    {
        "name": "A6_3layer_genebatch",
        "n_latent": 30, "n_hidden": 256, "n_layers": 3,
        "gene_likelihood": "nb", "dispersion": "gene-batch",
        "classification_ratio": 100,
        "pretrain_epochs": 100, "finetune_epochs": 50,
        "label_key": "celltype_broad",
        "n_samples_per_label": None,
        "plan_kwargs_extra": {},
    },
    {
        "name": "A7_lat50",
        "n_latent": 50, "n_hidden": 256, "n_layers": 2,
        "gene_likelihood": "nb", "dispersion": "gene",
        "classification_ratio": 100,
        "pretrain_epochs": 100, "finetune_epochs": 50,
        "label_key": "celltype_broad",
        "n_samples_per_label": None,
        "plan_kwargs_extra": {},
    },
    {
        "name": "A8_subsample",
        "n_latent": 30, "n_hidden": 256, "n_layers": 2,
        "gene_likelihood": "nb", "dispersion": "gene",
        "classification_ratio": 200,
        "pretrain_epochs": 100, "finetune_epochs": 50,
        "label_key": "celltype_broad",
        "n_samples_per_label": 100,
        "plan_kwargs_extra": {},
    },
    {
        "name": "A9_fine_small",
        "n_latent": 20, "n_hidden": 128, "n_layers": 2,
        "gene_likelihood": "zinb", "dispersion": "gene",
        "classification_ratio": 100,
        "pretrain_epochs": 50, "finetune_epochs": 30,
        "label_key": "celltype",
        "n_samples_per_label": None,
        "plan_kwargs_extra": {},
    },
    {
        "name": "A10_fine_big",
        "n_latent": 30, "n_hidden": 256, "n_layers": 2,
        "gene_likelihood": "nb", "dispersion": "gene",
        "classification_ratio": 100,
        "pretrain_epochs": 100, "finetune_epochs": 50,
        "label_key": "celltype",
        "n_samples_per_label": None,
        "plan_kwargs_extra": {},
    },
    {
        "name": "A11_kl_warmup",
        "n_latent": 30, "n_hidden": 256, "n_layers": 2,
        "gene_likelihood": "nb", "dispersion": "gene",
        "classification_ratio": 100,
        "pretrain_epochs": 100, "finetune_epochs": 50,
        "label_key": "celltype_broad",
        "n_samples_per_label": None,
        "plan_kwargs_extra": {"n_epochs_kl_warmup": 15},
    },
    {
        "name": "A12_long_train",
        "n_latent": 30, "n_hidden": 256, "n_layers": 2,
        "gene_likelihood": "nb", "dispersion": "gene",
        "classification_ratio": 50,
        "pretrain_epochs": 200, "finetune_epochs": 80,
        "label_key": "celltype_broad",
        "n_samples_per_label": None,
        "plan_kwargs_extra": {},
    },
]


# ============================================================
# Metric helpers
# ============================================================

def compute_metrics(emb, obs, batch_key="platform"):
    """Compute batch_asw, batch_score, ct_broad_asw, ct_asw, platform_entropy."""
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

    # Build a temporary AnnData just for neighbors/leiden
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
    print(f"  scANVI-max sweep: {len(CONFIGS)} configurations")
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
    if "celltype" in adata_orig.obs:
        n_ct = adata_orig.obs["celltype"].nunique()
        print(f"  celltype: {n_ct} unique types")

    results = []
    best_ct_broad_asw = -1.0
    best_config_name = None
    best_embedding = None

    for i, cfg in enumerate(CONFIGS):
        name = cfg["name"]
        print(f"\n{'='*70}")
        print(f"  [{i+1}/{len(CONFIGS)}] {name}")
        print(f"  n_latent={cfg['n_latent']}, n_hidden={cfg['n_hidden']}, "
              f"n_layers={cfg['n_layers']}, gene_likelihood={cfg['gene_likelihood']}, "
              f"dispersion={cfg['dispersion']}, cls_ratio={cfg['classification_ratio']}, "
              f"pretrain={cfg['pretrain_epochs']}, finetune={cfg['finetune_epochs']}, "
              f"label_key={cfg['label_key']}, n_samples_per_label={cfg['n_samples_per_label']}")
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
                adata, n_top_genes=2000, batch_key="platform", flavor="seurat_v3"
            )
            adata = adata[:, adata.var["highly_variable"]].copy()
            print(f"  HVG subset: {adata.n_vars} genes")

            # Prepare labels
            label_key = cfg["label_key"]
            adata.obs["_scanvi_labels"] = (
                adata.obs[label_key].astype(str).replace("nan", "Unknown")
            )
            n_labeled = (adata.obs["_scanvi_labels"] != "Unknown").sum()
            n_types = adata.obs["_scanvi_labels"].nunique()
            print(f"  Labels ({label_key}): {n_labeled}/{adata.n_obs} labeled, {n_types} types")

            # --- scVI pretrain ---
            print(f"  Training scVI pretrain ({cfg['pretrain_epochs']} epochs)...")
            scvi.model.SCVI.setup_anndata(adata, batch_key="platform")
            scvi_model = scvi.model.SCVI(
                adata,
                n_latent=cfg["n_latent"],
                n_hidden=cfg["n_hidden"],
                n_layers=cfg["n_layers"],
                gene_likelihood=cfg["gene_likelihood"],
                dispersion=cfg["dispersion"],
            )

            t0 = time.time()
            scvi_model.train(
                max_epochs=cfg["pretrain_epochs"],
                early_stopping=True,
                early_stopping_patience=10,
                batch_size=512,
                train_size=0.9,
            )
            pretrain_time = time.time() - t0
            print(f"  scVI pretrain done ({pretrain_time:.0f}s)")

            # --- scANVI finetune ---
            print(f"  Training scANVI ({cfg['finetune_epochs']} epochs, cls_ratio={cfg['classification_ratio']})...")

            plan_kwargs = {"classification_ratio": cfg["classification_ratio"]}
            plan_kwargs.update(cfg["plan_kwargs_extra"])

            scanvi_kwargs = {
                "unlabeled_category": "Unknown",
                "labels_key": "_scanvi_labels",
            }
            if cfg["n_samples_per_label"] is not None:
                scanvi_kwargs["n_samples_per_label"] = cfg["n_samples_per_label"]

            scanvi_model = scvi.model.SCANVI.from_scvi_model(
                scvi_model, **scanvi_kwargs,
            )

            t0 = time.time()
            scanvi_model.train(
                max_epochs=cfg["finetune_epochs"],
                batch_size=512,
                train_size=0.9,
                plan_kwargs=plan_kwargs,
            )
            finetune_time = time.time() - t0
            print(f"  scANVI finetune done ({finetune_time:.0f}s)")

            # Get latent + predictions
            latent = scanvi_model.get_latent_representation()
            predictions = scanvi_model.predict()
            print(f"  Latent: {latent.shape}")

            # Prediction accuracy
            orig = adata.obs[label_key].astype(str)
            pred = pd.Series(predictions, index=adata.obs.index).astype(str)
            valid = (orig != "nan") & (orig != "Unknown")
            acc = (orig[valid] == pred[valid]).mean() if valid.sum() > 0 else None
            if acc is not None:
                print(f"  Prediction accuracy: {acc:.3f} ({valid.sum()} cells)")

            # Compute metrics
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
                "pred_accuracy": round(acc, 4) if acc is not None else None,
                "pretrain_time_s": round(pretrain_time, 1),
                "finetune_time_s": round(finetune_time, 1),
                "total_time_s": round(total_time, 1),
            }
            results.append(row)

            print(f"\n  RESULTS: batch_score={batch_score:.4f}, "
                  f"ct_broad_asw={ct_broad_asw:.4f if ct_broad_asw is not None else 'N/A'}, "
                  f"ct_asw={ct_asw:.4f if ct_asw is not None else 'N/A'}, "
                  f"entropy={platform_entropy:.4f}, "
                  f"acc={acc:.4f if acc is not None else 'N/A'}")

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
                "pretrain_time_s": None, "finetune_time_s": None,
                "total_time_s": None,
            })
            continue

    # ============================================================
    # Final summary
    # ============================================================
    print(f"\n{'='*70}")
    print("  SCANVI-MAX SWEEP COMPLETE")
    print(f"{'='*70}")

    if results:
        results_df = pd.DataFrame(results)
        results_df_sorted = results_df.sort_values("ct_broad_asw", ascending=False)

        print("\nFinal ranking (by ct_broad_asw descending):")
        print(results_df_sorted.to_string(index=False))

        csv_path = OUTDIR / "scanvi_max_sweep.csv"
        results_df_sorted.to_csv(csv_path, index=False)
        print(f"\nSaved: {csv_path}")

    if best_embedding is not None:
        print(f"\nBest config: {best_config_name} (ct_broad_asw={best_ct_broad_asw:.4f})")
        print(f"Saving adata with X_scanvi_max embedding...")

        adata_out = adata_orig.copy()
        adata_out.obsm["X_scanvi_max"] = best_embedding

        out_path = OUTDIR / "adata.m2.scanvi_max.h5ad"
        adata_out.write_h5ad(out_path)
        print(f"Saved: {out_path} ({adata_out.n_obs} cells)")
    else:
        print("\nNo successful configs — no adata saved.")

    print("\nDone.")


if __name__ == "__main__":
    main()
