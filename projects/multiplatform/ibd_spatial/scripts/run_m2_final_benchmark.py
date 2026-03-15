#!/usr/bin/env python3
"""M2 final benchmark: run top scANVI configs, compute UMAPs, generate report.

Loads adata.m2.h5ad (which already has embeddings for PCA, Harmony, scVI,
Scanorama, scANVI-vanilla, resolVI, resolVI-SS) and adds optimized scANVI
configs from the hyperparameter search.

For each method, stores:
  - obsm['X_{method}'] — latent embedding
  - obsm['X_umap_{method}'] — 2D UMAP

Saves per-method lightweight files (obs + obsm, no X) to
results/m2_benchmark/embeddings/{method}.h5ad for downstream report generation.

Generates:
  - results/m2_benchmark/m2_final_benchmark.csv — full metrics table
  - figures/m2/m2_umap_grid.png — summary UMAP grid
  - figures/m2/m2_umap_{method}.png — per-method UMAP quartet
  - figures/QC/m2_final_benchmark_report.html — standalone HTML report

Usage:
    python run_m2_final_benchmark.py
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
from matplotlib.gridspec import GridSpec
from sklearn.metrics import silhouette_score
from scipy.stats import entropy as scipy_entropy

import torch

torch.set_float32_matmul_precision("medium")

import scvi

WORKDIR = Path("/home/fs01/juk4007/elementolab/projects/ibd_spatial")
INPUT_PATH = WORKDIR / "results" / "m2_benchmark" / "adata.m2.h5ad"
OUTDIR = WORKDIR / "results" / "m2_benchmark"
EMBDIR = OUTDIR / "embeddings"
FIGDIR = WORKDIR / "figures" / "m2"
REPORTDIR = WORKDIR / "figures" / "QC"

EMBDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)
REPORTDIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# Top scANVI configs from random search (2026-03-13)
# ============================================================
# Fixed params (confirmed optimal)
FIXED = {
    "gene_likelihood": "nb",
    "dispersion": "gene-batch",
    "label_key": "celltype_broad",
    "batch_size": 512,
    "train_size": 0.9,
}

SCANVI_CONFIGS = {
    "scanvi_r023": {
        "n_layers": 4, "n_latent": 53, "n_hidden": 512,
        "classification_ratio": 192, "dropout_rate": 0.20,
        "pretrain_epochs": 125, "finetune_epochs": 30,
        "description": "Random search winner (ct_broad_asw=0.396)",
    },
    "scanvi_r026": {
        "n_layers": 4, "n_latent": 42, "n_hidden": 128,
        "classification_ratio": 178, "dropout_rate": 0.10,
        "pretrain_epochs": 75, "finetune_epochs": 75,
        "description": "Runner-up (ct_broad_asw=0.374)",
    },
    "scanvi_a6": {
        "n_layers": 3, "n_latent": 30, "n_hidden": 256,
        "classification_ratio": 100, "dropout_rate": 0.10,
        "pretrain_epochs": 100, "finetune_epochs": 50,
        "description": "Grid sweep best (ct_broad_asw=0.189)",
    },
}

# Existing methods already in adata.m2.h5ad (just need UMAPs if missing)
EXISTING_METHODS = ["pca", "harmony", "scvi", "scanorama", "scanvi", "resolvi", "resolvi_ss"]

# Color keys for UMAP grid
COLOR_KEYS = ["platform", "celltype_broad", "disease", "patient_id"]

# ============================================================
# Metric computation
# ============================================================

def compute_all_metrics(emb, obs, batch_key="platform"):
    """Compute batch and bio metrics for one embedding."""
    n_sample = min(10_000, emb.shape[0])

    # Batch ASW
    batch_asw = silhouette_score(
        emb, obs[batch_key].values,
        sample_size=n_sample, random_state=42,
    )
    batch_score = 1 - abs(batch_asw)

    # Celltype broad ASW
    ct_broad_asw = None
    labels = obs["celltype_broad"].values
    valid = pd.notna(labels)
    if valid.sum() > 100 and len(np.unique(labels[valid])) > 1:
        ct_broad_asw = silhouette_score(
            emb[valid], labels[valid],
            sample_size=min(n_sample, int(valid.sum())), random_state=42,
        )

    # Celltype ASW
    ct_asw = None
    if "celltype" in obs.columns:
        labels_ct = obs["celltype"].values
        valid_ct = pd.notna(labels_ct)
        if valid_ct.sum() > 100 and len(np.unique(labels_ct[valid_ct])) > 1:
            ct_asw = silhouette_score(
                emb[valid_ct], labels_ct[valid_ct],
                sample_size=min(n_sample, int(valid_ct.sum())), random_state=42,
            )

    # Platform entropy
    import warnings
    warnings.filterwarnings("ignore", category=FutureWarning)
    tmp = ad.AnnData(obs=obs.copy())
    tmp.obsm["X_emb"] = emb
    sc.pp.neighbors(tmp, use_rep="X_emb")
    sc.tl.leiden(tmp, resolution=0.5, key_added="leiden_tmp")
    n_batches = obs[batch_key].nunique()
    max_ent = np.log(n_batches) if n_batches > 1 else 1.0
    entropies = []
    for cluster in tmp.obs["leiden_tmp"].unique():
        mask = tmp.obs["leiden_tmp"] == cluster
        counts = obs.loc[mask, batch_key].value_counts()
        ent = scipy_entropy(counts / counts.sum()) / max_ent if max_ent > 0 else 0
        entropies.append(ent)
    platform_entropy = np.median(entropies)

    return {
        "batch_asw": round(batch_asw, 4),
        "batch_score": round(batch_score, 4),
        "ct_broad_asw": round(ct_broad_asw, 4) if ct_broad_asw is not None else None,
        "ct_asw": round(ct_asw, 4) if ct_asw is not None else None,
        "platform_entropy": round(platform_entropy, 4),
    }


# ============================================================
# scANVI training
# ============================================================

def run_scanvi(adata_full, config_name, params):
    """Train scVI + scANVI and return latent embedding."""
    print(f"\n{'=' * 70}")
    print(f"  Training {config_name}: {params['description']}")
    print(f"  n_layers={params['n_layers']}, n_latent={params['n_latent']}, "
          f"n_hidden={params['n_hidden']}, cls_ratio={params['classification_ratio']}, "
          f"dropout={params['dropout_rate']}")
    print(f"{'=' * 70}")

    adata = adata_full.copy()
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()

    # HVG
    sc.pp.highly_variable_genes(
        adata, n_top_genes=2000, batch_key="platform", flavor="seurat_v3"
    )
    adata = adata[:, adata.var["highly_variable"]].copy()
    print(f"  HVG subset: {adata.n_vars} genes")

    # Labels
    label_key = FIXED["label_key"]
    adata.obs["_scanvi_labels"] = (
        adata.obs[label_key].astype(str).replace("nan", "Unknown")
    )
    n_labeled = (adata.obs["_scanvi_labels"] != "Unknown").sum()
    print(f"  Labels: {n_labeled}/{adata.n_obs} labeled")

    # scVI pretrain
    print(f"  scVI pretrain ({params['pretrain_epochs']} epochs max)...")
    scvi.model.SCVI.setup_anndata(adata, batch_key="platform")
    scvi_model = scvi.model.SCVI(
        adata,
        n_latent=params["n_latent"],
        n_hidden=params["n_hidden"],
        n_layers=params["n_layers"],
        gene_likelihood=FIXED["gene_likelihood"],
        dispersion=FIXED["dispersion"],
        dropout_rate=params["dropout_rate"],
    )
    t0 = time.time()
    scvi_model.train(
        max_epochs=params["pretrain_epochs"],
        early_stopping=True,
        early_stopping_patience=10,
        batch_size=FIXED["batch_size"],
        train_size=FIXED["train_size"],
    )
    pretrain_time = time.time() - t0
    pretrain_actual = scvi_model.history["elbo_train"].shape[0]
    print(f"  Pretrain done: {pretrain_time:.0f}s, {pretrain_actual} epochs")

    # scANVI finetune
    print(f"  scANVI finetune ({params['finetune_epochs']} epochs)...")
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        unlabeled_category="Unknown",
        labels_key="_scanvi_labels",
    )
    t0 = time.time()
    scanvi_model.train(
        max_epochs=params["finetune_epochs"],
        batch_size=FIXED["batch_size"],
        train_size=FIXED["train_size"],
        plan_kwargs={"classification_ratio": params["classification_ratio"]},
    )
    finetune_time = time.time() - t0
    finetune_actual = scanvi_model.history["elbo_train"].shape[0]
    print(f"  Finetune done: {finetune_time:.0f}s, {finetune_actual} epochs")

    latent = scanvi_model.get_latent_representation()
    predictions = scanvi_model.predict()
    print(f"  Latent shape: {latent.shape}")

    # Prediction accuracy
    orig = adata_full.obs[label_key].astype(str)
    pred = pd.Series(predictions, index=adata.obs.index).astype(str)
    valid = (orig != "nan") & (orig != "Unknown")
    acc = (orig[valid] == pred[valid]).mean() if valid.sum() > 0 else None
    if acc is not None:
        print(f"  Prediction accuracy: {acc:.4f}")

    return latent, pretrain_time + finetune_time, acc


# ============================================================
# UMAP computation
# ============================================================

def compute_umap(adata, emb_key, umap_key, n_neighbors=15):
    """Compute UMAP from embedding and store in obsm."""
    if umap_key in adata.obsm:
        print(f"  UMAP {umap_key} already exists, skipping")
        return
    print(f"  Computing UMAP: {emb_key} -> {umap_key}...")
    # Use a temporary AnnData to avoid disturbing the main graph
    tmp = ad.AnnData(obs=adata.obs.copy())
    tmp.obsm["X_emb"] = adata.obsm[emb_key]
    sc.pp.neighbors(tmp, use_rep="X_emb", n_neighbors=n_neighbors)
    sc.tl.umap(tmp)
    adata.obsm[umap_key] = tmp.obsm["X_umap"].copy()
    print(f"  Done: {umap_key}")


# ============================================================
# Plotting
# ============================================================

def plot_umap_grid(adata, methods, color_keys, figdir):
    """Plot a grid of UMAPs: rows=methods, cols=color_keys."""
    n_rows = len(methods)
    n_cols = len(color_keys)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4.5 * n_rows))
    if n_rows == 1:
        axes = axes[np.newaxis, :]
    if n_cols == 1:
        axes = axes[:, np.newaxis]

    for i, method in enumerate(methods):
        umap_key = f"X_umap_{method}"
        if umap_key not in adata.obsm:
            for j in range(n_cols):
                axes[i, j].text(0.5, 0.5, f"No UMAP\n{method}",
                               ha="center", va="center", fontsize=12)
                axes[i, j].set_xlim(0, 1)
                axes[i, j].set_ylim(0, 1)
            continue

        for j, color_key in enumerate(color_keys):
            ax = axes[i, j]
            coords = adata.obsm[umap_key]

            # Subsample for speed
            n = min(50_000, coords.shape[0])
            idx = np.random.RandomState(42).choice(coords.shape[0], n, replace=False)

            categories = adata.obs[color_key].values[idx]
            unique_cats = sorted(set(str(c) for c in categories if pd.notna(c)))

            if len(unique_cats) <= 20:
                cmap = plt.cm.tab20 if len(unique_cats) > 10 else plt.cm.tab10
                color_map = {cat: cmap(k / max(len(unique_cats), 1))
                             for k, cat in enumerate(unique_cats)}
                colors = [color_map.get(str(c), (0.8, 0.8, 0.8, 1.0))
                          for c in categories]
            else:
                colors = "grey"

            ax.scatter(coords[idx, 0], coords[idx, 1], c=colors,
                       s=0.5, alpha=0.5, rasterized=True)
            ax.set_title(f"{method} | {color_key}", fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)

            # Legend for small number of categories
            if len(unique_cats) <= 10 and isinstance(colors, list):
                from matplotlib.lines import Line2D
                handles = [Line2D([0], [0], marker="o", color="w",
                                  markerfacecolor=color_map[cat], markersize=6,
                                  label=cat)
                           for cat in unique_cats]
                ax.legend(handles=handles, loc="upper right", fontsize=5,
                          framealpha=0.7, markerscale=0.8)

    plt.tight_layout()
    out_path = figdir / "m2_umap_grid.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")
    return out_path


def plot_single_method_umap(adata, method, color_keys, figdir):
    """Plot 4-panel UMAP for one method."""
    umap_key = f"X_umap_{method}"
    if umap_key not in adata.obsm:
        return None

    n_cols = len(color_keys)
    fig, axes = plt.subplots(1, n_cols, figsize=(5 * n_cols, 4.5))
    if n_cols == 1:
        axes = [axes]

    coords = adata.obsm[umap_key]
    n = min(50_000, coords.shape[0])
    idx = np.random.RandomState(42).choice(coords.shape[0], n, replace=False)

    for j, color_key in enumerate(color_keys):
        ax = axes[j]
        categories = adata.obs[color_key].values[idx]
        unique_cats = sorted(set(str(c) for c in categories if pd.notna(c)))

        cmap = plt.cm.tab20 if len(unique_cats) > 10 else plt.cm.tab10
        color_map = {cat: cmap(k / max(len(unique_cats), 1))
                     for k, cat in enumerate(unique_cats)}
        colors = [color_map.get(str(c), (0.8, 0.8, 0.8, 1.0))
                  for c in categories]

        ax.scatter(coords[idx, 0], coords[idx, 1], c=colors,
                   s=1.0, alpha=0.5, rasterized=True)
        ax.set_title(f"{color_key}", fontsize=11)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

        if len(unique_cats) <= 15:
            from matplotlib.lines import Line2D
            handles = [Line2D([0], [0], marker="o", color="w",
                              markerfacecolor=color_map[cat], markersize=6,
                              label=cat)
                       for cat in unique_cats]
            ax.legend(handles=handles, loc="upper right", fontsize=6,
                      framealpha=0.7, markerscale=0.8)

    fig.suptitle(method, fontsize=14, fontweight="bold")
    plt.tight_layout()
    out_path = figdir / f"m2_umap_{method}.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")
    return out_path


# ============================================================
# HTML report generation
# ============================================================

def generate_html_report(metrics_df, umap_grid_path, per_method_paths, report_path):
    """Generate standalone HTML benchmark report."""
    import base64

    def img_to_b64(path):
        if path and Path(path).exists():
            data = Path(path).read_bytes()
            return base64.b64encode(data).decode()
        return None

    # Sort by ct_broad_asw descending
    df = metrics_df.sort_values("ct_broad_asw", ascending=False).reset_index(drop=True)

    # Build metrics table HTML
    table_rows = []
    for _, row in df.iterrows():
        is_best = row.name == 0
        style = ' style="background:#e8f5e9;font-weight:bold;"' if is_best else ""
        table_rows.append(f"""<tr{style}>
            <td>{row.name + 1}</td>
            <td>{row['method']}</td>
            <td>{row['ct_broad_asw']:.4f}</td>
            <td>{row['batch_score']:.4f}</td>
            <td>{row.get('platform_entropy', 'N/A')}</td>
            <td>{row.get('ct_asw', 'N/A')}</td>
            <td>{row.get('pred_accuracy', 'N/A')}</td>
            <td>{row.get('description', '')}</td>
        </tr>""")

    table_html = "\n".join(table_rows)

    # Embed images
    grid_b64 = img_to_b64(umap_grid_path)
    grid_img = f'<img src="data:image/png;base64,{grid_b64}" style="max-width:100%;">' if grid_b64 else "<p>No grid image</p>"

    method_imgs = []
    for method, path in per_method_paths.items():
        b64 = img_to_b64(path)
        if b64:
            method_imgs.append(f"""
            <div class="method-panel">
                <h3>{method}</h3>
                <img src="data:image/png;base64,{b64}" style="max-width:100%;">
            </div>""")

    method_imgs_html = "\n".join(method_imgs)

    # Best method summary
    best = df.iloc[0]

    html = f"""<!DOCTYPE html>
<html><head>
<meta charset="utf-8">
<title>M2 Integration Benchmark Report</title>
<style>
    body {{ font-family: Helvetica, Arial, sans-serif; max-width: 1400px; margin: 0 auto; padding: 20px; background: #fafafa; }}
    h1 {{ color: #1a237e; border-bottom: 3px solid #1a237e; padding-bottom: 10px; }}
    h2 {{ color: #283593; margin-top: 30px; }}
    .summary-card {{ background: #e8eaf6; border-radius: 8px; padding: 15px 20px; margin: 15px 0; display: inline-block; min-width: 200px; }}
    .summary-card .value {{ font-size: 2em; font-weight: bold; color: #1a237e; }}
    .summary-card .label {{ font-size: 0.9em; color: #5c6bc0; }}
    table {{ border-collapse: collapse; width: 100%; margin: 15px 0; }}
    th {{ background: #283593; color: white; padding: 10px 12px; text-align: left; font-size: 0.9em; }}
    td {{ padding: 8px 12px; border-bottom: 1px solid #e0e0e0; font-size: 0.9em; }}
    tr:hover {{ background: #f5f5f5; }}
    .method-panel {{ margin: 20px 0; padding: 15px; background: white; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}
    .grid-container {{ background: white; padding: 20px; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); margin: 20px 0; }}
    .footer {{ color: #9e9e9e; font-size: 0.8em; margin-top: 40px; border-top: 1px solid #e0e0e0; padding-top: 10px; }}
</style>
</head><body>

<h1>M2 Cross-Platform Integration Benchmark</h1>
<p><strong>Dataset:</strong> CosMx 6k + Xenium 5K | 8 samples, 4 patients | 164,392 cells | 2,552 shared genes</p>
<p><strong>Primary metric:</strong> ct_broad_asw (bio conservation via celltype_broad silhouette)</p>

<div>
    <div class="summary-card">
        <div class="value">{best['ct_broad_asw']:.4f}</div>
        <div class="label">Best ct_broad ASW ({best['method']})</div>
    </div>
    <div class="summary-card">
        <div class="value">{best['batch_score']:.4f}</div>
        <div class="label">Batch Score ({best['method']})</div>
    </div>
    <div class="summary-card">
        <div class="value">{len(df)}</div>
        <div class="label">Methods Compared</div>
    </div>
</div>

<h2>Integration Comparison</h2>
<table>
    <tr>
        <th>Rank</th><th>Method</th><th>ct_broad ASW</th><th>Batch Score</th>
        <th>Platform Entropy</th><th>ct ASW</th><th>Pred Accuracy</th><th>Notes</th>
    </tr>
    {table_html}
</table>

<h2>UMAP Grid (all methods)</h2>
<div class="grid-container">
    {grid_img}
</div>

<h2>Per-Method UMAPs</h2>
{method_imgs_html}

<div class="footer">
    Generated by sc_tools M2 benchmark pipeline.
    Dataset: IBD spatial (CosMx 6k + Xenium 5K, patient-matched).
</div>
</body></html>"""

    Path(report_path).write_text(html)
    print(f"Report saved: {report_path}")


# ============================================================
# Main
# ============================================================

def main():
    t_start = time.time()

    # --- 1. Load adata ---
    print("=== Loading adata.m2.h5ad ===")
    # Copy to /tmp for faster I/O
    local_path = Path(f"/tmp/adata_benchmark_{os.getpid()}.h5ad")
    if not local_path.exists():
        import shutil
        print(f"  Copying to {local_path}...")
        shutil.copy2(INPUT_PATH, local_path)
    adata = ad.read_h5ad(local_path)
    print(f"  Shape: {adata.n_obs} x {adata.n_vars}")
    print(f"  Existing obsm: {list(adata.obsm.keys())}")

    # --- 2. Train new scANVI configs ---
    print("\n=== Training optimized scANVI configs ===")
    training_info = {}
    for config_name, params in SCANVI_CONFIGS.items():
        emb_key = f"X_{config_name}"
        if emb_key in adata.obsm:
            print(f"\n  {config_name} already in obsm, skipping training")
            training_info[config_name] = {"time": 0, "pred_accuracy": None}
            continue

        latent, train_time, pred_acc = run_scanvi(adata, config_name, params)
        adata.obsm[emb_key] = latent
        training_info[config_name] = {"time": train_time, "pred_accuracy": pred_acc}
        print(f"  Stored: obsm['{emb_key}'] shape={latent.shape}")

    # --- 3. Compute UMAPs for all methods ---
    print("\n=== Computing UMAPs ===")
    all_methods = EXISTING_METHODS + list(SCANVI_CONFIGS.keys())
    for method in all_methods:
        emb_key = f"X_{method}"
        umap_key = f"X_umap_{method}"
        if emb_key in adata.obsm:
            compute_umap(adata, emb_key, umap_key)
        else:
            print(f"  SKIP {method}: no embedding in obsm")

    # --- 4. Compute metrics for all methods ---
    print("\n=== Computing metrics ===")
    all_results = []
    for method in all_methods:
        emb_key = f"X_{method}"
        if emb_key not in adata.obsm:
            continue

        print(f"\n  Metrics for {method}...")
        metrics = compute_all_metrics(adata.obsm[emb_key], adata.obs)
        info = training_info.get(method, {})

        row = {
            "method": method,
            **metrics,
            "pred_accuracy": info.get("pred_accuracy"),
            "n_latent_dims": adata.obsm[emb_key].shape[1],
        }

        # Add description for scANVI configs
        if method in SCANVI_CONFIGS:
            row["description"] = SCANVI_CONFIGS[method]["description"]
        elif method == "scanvi":
            row["description"] = "Vanilla scANVI (stock params)"
        elif method == "scvi":
            row["description"] = "scVI (unsupervised)"
        elif method == "harmony":
            row["description"] = "Harmony (linear correction)"
        elif method == "resolvi":
            row["description"] = "resolVI (spatial-aware unsupervised)"
        elif method == "resolvi_ss":
            row["description"] = "resolVI semi-supervised"
        else:
            row["description"] = ""

        all_results.append(row)
        print(f"    batch_score={metrics['batch_score']}, "
              f"ct_broad_asw={metrics['ct_broad_asw']}")

    results_df = pd.DataFrame(all_results)
    csv_path = OUTDIR / "m2_final_benchmark.csv"
    results_df.to_csv(csv_path, index=False)
    print(f"\nSaved metrics: {csv_path}")
    print(results_df.sort_values("ct_broad_asw", ascending=False).to_string(index=False))

    # --- 5. Save per-method lightweight embeddings ---
    print("\n=== Saving per-method embedding files ===")
    for method in all_methods:
        emb_key = f"X_{method}"
        umap_key = f"X_umap_{method}"
        if emb_key not in adata.obsm:
            continue

        # Lightweight: obs + obsm only (no X, no var)
        emb_adata = ad.AnnData(obs=adata.obs.copy())
        emb_adata.obsm[emb_key] = adata.obsm[emb_key]
        if umap_key in adata.obsm:
            emb_adata.obsm[umap_key] = adata.obsm[umap_key]
        emb_path = EMBDIR / f"{method}.h5ad"
        emb_adata.write_h5ad(emb_path)
        size_mb = emb_path.stat().st_size / 1e6
        print(f"  {method}: {emb_path} ({size_mb:.1f} MB)")

    # --- 6. Plot UMAPs ---
    print("\n=== Plotting UMAPs ===")
    methods_with_umap = [m for m in all_methods if f"X_umap_{m}" in adata.obsm]
    grid_path = plot_umap_grid(adata, methods_with_umap, COLOR_KEYS, FIGDIR)

    per_method_paths = {}
    for method in methods_with_umap:
        path = plot_single_method_umap(adata, method, COLOR_KEYS, FIGDIR)
        if path:
            per_method_paths[method] = path

    # --- 7. Generate HTML report ---
    print("\n=== Generating HTML report ===")
    report_path = REPORTDIR / "m2_final_benchmark_report.html"
    generate_html_report(results_df, grid_path, per_method_paths, report_path)

    # --- 8. Save updated adata with all embeddings ---
    print("\n=== Saving updated adata ===")
    out_path = OUTDIR / "adata.m2.h5ad"
    adata.write_h5ad(out_path)
    size_gb = out_path.stat().st_size / 1e9
    print(f"Saved: {out_path} ({size_gb:.2f} GB)")

    # Cleanup temp file
    if local_path.exists():
        local_path.unlink()

    total_time = time.time() - t_start
    print(f"\n{'=' * 70}")
    print(f"  COMPLETE: {total_time:.0f}s ({total_time / 60:.1f} min)")
    print(f"  Methods benchmarked: {len(results_df)}")
    print(f"  Best: {results_df.sort_values('ct_broad_asw', ascending=False).iloc[0]['method']} "
          f"(ct_broad_asw={results_df['ct_broad_asw'].max():.4f})")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
