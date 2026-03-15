#!/usr/bin/env python3
"""scANVI random hyperparameter search: parallelized via SLURM job arrays.

Each task runs ONE scANVI config (scVI pretrain + scANVI finetune) and saves
metrics to a per-task CSV. The SLURM_ARRAY_TASK_ID maps to a config index.

Design:
- 40 configs total (R000-R039), sparse random search over a broader space
- Config R000 = A6 baseline (n_layers=3, n_latent=30, n_hidden=256,
  classification_ratio=100, pretrain=100, finetune=50, dropout=0.1)
- Fixed (confirmed optimal): nb, gene-batch, celltype_broad, batch_size=512
- Random search params: n_layers, n_latent, n_hidden, classification_ratio,
  pretrain_epochs, finetune_epochs, dropout_rate
- Fixed seed (42) for reproducibility

Reads:  results/m2_benchmark/adata.m2.h5ad (164K cells, 2552 genes)
Writes: results/scanvi_rsearch/config_{idx:03d}.csv
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
OUTDIR = WORKDIR / "results" / "scanvi_rsearch"

# ============================================================
# Fixed parameters (confirmed optimal from A1-A12 sweep)
# ============================================================
FIXED = {
    "gene_likelihood": "nb",
    "dispersion": "gene-batch",
    "label_key": "celltype_broad",
    "batch_size": 512,
    "train_size": 0.9,
}

# ============================================================
# Random search space
# ============================================================
SEARCH_SPACE = {
    "n_layers": [2, 3, 4],
    "n_latent": (15, 60),           # uniform_int
    "n_hidden": [128, 192, 256, 384, 512],
    "classification_ratio": (50, 200),  # uniform_int
    "pretrain_epochs": [75, 100, 125, 150],
    "finetune_epochs": [30, 50, 75],
    "dropout_rate": [0.05, 0.1, 0.15, 0.2],
}

# A6 baseline config
A6_BASELINE = {
    "n_layers": 3,
    "n_latent": 30,
    "n_hidden": 256,
    "classification_ratio": 100,
    "pretrain_epochs": 100,
    "finetune_epochs": 50,
    "dropout_rate": 0.1,
}

N_CONFIGS = 40
SEED = 42


def build_config_list(n_configs=N_CONFIGS, seed=SEED):
    """Build a random search of configs from the search space.

    Config R000 is always the A6 baseline. R001-R039 are random samples
    drawn with a fixed seed for reproducibility.
    """
    rng = np.random.RandomState(seed)

    configs = []

    # R000: A6 baseline
    cfg0 = dict(A6_BASELINE)
    cfg0.update(FIXED)
    cfg0["config_idx"] = 0
    cfg0["name"] = (
        f"R000_L{cfg0['n_layers']}_d{cfg0['n_latent']}_h{cfg0['n_hidden']}"
        f"_c{cfg0['classification_ratio']}_p{cfg0['pretrain_epochs']}"
        f"_f{cfg0['finetune_epochs']}_dr{cfg0['dropout_rate']}"
    )
    configs.append(cfg0)

    # R001-R039: random samples
    seen = {(A6_BASELINE["n_layers"], A6_BASELINE["n_latent"],
             A6_BASELINE["n_hidden"], A6_BASELINE["classification_ratio"],
             A6_BASELINE["pretrain_epochs"], A6_BASELINE["finetune_epochs"],
             A6_BASELINE["dropout_rate"])}

    attempts = 0
    while len(configs) < n_configs and attempts < 1000:
        attempts += 1

        n_layers = rng.choice(SEARCH_SPACE["n_layers"])
        n_latent = rng.randint(SEARCH_SPACE["n_latent"][0],
                               SEARCH_SPACE["n_latent"][1] + 1)
        n_hidden = rng.choice(SEARCH_SPACE["n_hidden"])
        classification_ratio = rng.randint(SEARCH_SPACE["classification_ratio"][0],
                                           SEARCH_SPACE["classification_ratio"][1] + 1)
        pretrain_epochs = rng.choice(SEARCH_SPACE["pretrain_epochs"])
        finetune_epochs = rng.choice(SEARCH_SPACE["finetune_epochs"])
        dropout_rate = rng.choice(SEARCH_SPACE["dropout_rate"])

        key = (n_layers, n_latent, n_hidden, classification_ratio,
               pretrain_epochs, finetune_epochs, dropout_rate)
        if key in seen:
            continue
        seen.add(key)

        idx = len(configs)
        cfg = {
            "n_layers": int(n_layers),
            "n_latent": int(n_latent),
            "n_hidden": int(n_hidden),
            "classification_ratio": int(classification_ratio),
            "pretrain_epochs": int(pretrain_epochs),
            "finetune_epochs": int(finetune_epochs),
            "dropout_rate": float(dropout_rate),
        }
        cfg.update(FIXED)
        cfg["config_idx"] = idx
        cfg["name"] = (
            f"R{idx:03d}_L{cfg['n_layers']}_d{cfg['n_latent']}_h{cfg['n_hidden']}"
            f"_c{cfg['classification_ratio']}_p{cfg['pretrain_epochs']}"
            f"_f{cfg['finetune_epochs']}_dr{cfg['dropout_rate']}"
        )
        configs.append(cfg)

    print(f"Generated {len(configs)} configs (seed={seed})")
    return configs


ALL_CONFIGS = build_config_list(n_configs=N_CONFIGS, seed=SEED)


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
# Run one config
# ============================================================
def run_one_config(config_idx: int):
    """Run a single scANVI config and save results."""
    cfg = ALL_CONFIGS[config_idx]
    name = cfg["name"]

    OUTDIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUTDIR / f"config_{config_idx:03d}.csv"

    # Check if already done
    if csv_path.exists():
        print(f"SKIP: {csv_path} already exists")
        return

    print(f"{'=' * 80}")
    print(f"  Config {config_idx}/{N_CONFIGS}: {name}")
    print(f"  n_layers={cfg['n_layers']}, n_latent={cfg['n_latent']}, "
          f"n_hidden={cfg['n_hidden']}, cls_ratio={cfg['classification_ratio']}, "
          f"pretrain={cfg['pretrain_epochs']}, finetune={cfg['finetune_epochs']}, "
          f"dropout={cfg['dropout_rate']}")
    print(f"  Fixed: likelihood={cfg['gene_likelihood']}, "
          f"dispersion={cfg['dispersion']}, labels={cfg['label_key']}")
    if config_idx == 0:
        print(f"  *** A6 BASELINE CONFIG ***")
    print(f"{'=' * 80}")

    # Load data (prefer override path, then /tmp copy, then original)
    override_path = os.environ.get("SCANVI_INPUT_OVERRIDE")
    local_path = Path("/tmp/adata.m2.h5ad")
    if override_path and Path(override_path).exists():
        load_path = Path(override_path)
    elif local_path.exists():
        load_path = local_path
    else:
        load_path = INPUT_PATH
    print(f"Loading adata from {load_path}...")
    adata_orig = ad.read_h5ad(load_path)
    print(f"  Shape: {adata_orig.n_obs} cells x {adata_orig.n_vars} genes")

    try:
        t0_total = time.time()

        # Prepare adata
        adata = adata_orig.copy()
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
        print(f"  Labels: {n_labeled}/{adata.n_obs} labeled")

        # --- scVI pretrain ---
        print(f"  scVI pretrain ({cfg['pretrain_epochs']} epochs max)...")
        scvi.model.SCVI.setup_anndata(adata, batch_key="platform")
        scvi_model = scvi.model.SCVI(
            adata,
            n_latent=cfg["n_latent"],
            n_hidden=cfg["n_hidden"],
            n_layers=cfg["n_layers"],
            gene_likelihood=cfg["gene_likelihood"],
            dispersion=cfg["dispersion"],
            dropout_rate=cfg["dropout_rate"],
        )

        t0 = time.time()
        scvi_model.train(
            max_epochs=cfg["pretrain_epochs"],
            early_stopping=True,
            early_stopping_patience=10,
            batch_size=cfg["batch_size"],
            train_size=cfg["train_size"],
        )
        pretrain_time = time.time() - t0
        pretrain_epochs_actual = scvi_model.history["elbo_train"].shape[0]
        print(f"  scVI pretrain done: {pretrain_time:.0f}s, "
              f"{pretrain_epochs_actual} epochs (of {cfg['pretrain_epochs']})")

        # --- scANVI finetune ---
        print(f"  scANVI finetune ({cfg['finetune_epochs']} epochs max, "
              f"cls_ratio={cfg['classification_ratio']})...")

        scanvi_model = scvi.model.SCANVI.from_scvi_model(
            scvi_model,
            unlabeled_category="Unknown",
            labels_key="_scanvi_labels",
        )

        plan_kwargs = {"classification_ratio": cfg["classification_ratio"]}

        t0 = time.time()
        scanvi_model.train(
            max_epochs=cfg["finetune_epochs"],
            batch_size=cfg["batch_size"],
            train_size=cfg["train_size"],
            plan_kwargs=plan_kwargs,
        )
        finetune_time = time.time() - t0
        finetune_epochs_actual = scanvi_model.history["elbo_train"].shape[0]
        print(f"  scANVI finetune done: {finetune_time:.0f}s, "
              f"{finetune_epochs_actual} epochs (of {cfg['finetune_epochs']})")

        # Get latent representation and predictions
        latent = scanvi_model.get_latent_representation()
        predictions = scanvi_model.predict()
        print(f"  Latent shape: {latent.shape}")

        # Prediction accuracy
        orig = adata.obs[label_key].astype(str)
        pred = pd.Series(predictions, index=adata.obs.index).astype(str)
        valid = (orig != "nan") & (orig != "Unknown")
        acc = (orig[valid] == pred[valid]).mean() if valid.sum() > 0 else None
        if acc is not None:
            print(f"  Prediction accuracy: {acc:.4f}")

        # Compute metrics
        print("  Computing metrics...")
        batch_asw, batch_score, ct_broad_asw, ct_asw = compute_metrics(
            latent, adata_orig.obs, batch_key="platform"
        )
        platform_entropy = compute_platform_entropy(
            latent, adata_orig.obs, batch_key="platform"
        )

        total_time = time.time() - t0_total

        # Build results row
        row = {
            "config_idx": config_idx,
            "name": name,
            "n_layers": cfg["n_layers"],
            "n_latent": cfg["n_latent"],
            "n_hidden": cfg["n_hidden"],
            "classification_ratio": cfg["classification_ratio"],
            "dropout_rate": cfg["dropout_rate"],
            "pretrain_epochs_max": cfg["pretrain_epochs"],
            "finetune_epochs_max": cfg["finetune_epochs"],
            "pretrain_epochs_actual": pretrain_epochs_actual,
            "finetune_epochs_actual": finetune_epochs_actual,
            "batch_asw": round(batch_asw, 4),
            "batch_score": round(batch_score, 4),
            "ct_broad_asw": round(ct_broad_asw, 4) if ct_broad_asw is not None else None,
            "ct_asw": round(ct_asw, 4) if ct_asw is not None else None,
            "platform_entropy": round(platform_entropy, 4),
            "pred_accuracy": round(acc, 4) if acc is not None else None,
            "pretrain_time_s": round(pretrain_time, 1),
            "finetune_time_s": round(finetune_time, 1),
            "total_time_s": round(total_time, 1),
            "status": "success",
        }

        print(f"\n  RESULTS:")
        print(f"    batch_score   = {batch_score:.4f}")
        print(f"    ct_broad_asw  = {ct_broad_asw:.4f}" if ct_broad_asw is not None else "    ct_broad_asw  = N/A")
        print(f"    ct_asw        = {ct_asw:.4f}" if ct_asw is not None else "    ct_asw        = N/A")
        print(f"    entropy       = {platform_entropy:.4f}")
        print(f"    pred_accuracy = {acc:.4f}" if acc is not None else "    pred_accuracy = N/A")
        print(f"    total_time    = {total_time:.0f}s")

    except Exception as e:
        print(f"  FAILED: {e}")
        traceback.print_exc()
        total_time = time.time() - t0_total if "t0_total" in dir() else 0
        row = {
            "config_idx": config_idx,
            "name": name,
            "n_layers": cfg["n_layers"],
            "n_latent": cfg["n_latent"],
            "n_hidden": cfg["n_hidden"],
            "classification_ratio": cfg["classification_ratio"],
            "dropout_rate": cfg["dropout_rate"],
            "pretrain_epochs_max": cfg["pretrain_epochs"],
            "finetune_epochs_max": cfg["finetune_epochs"],
            "pretrain_epochs_actual": None,
            "finetune_epochs_actual": None,
            "batch_asw": None,
            "batch_score": None,
            "ct_broad_asw": None,
            "ct_asw": None,
            "platform_entropy": None,
            "pred_accuracy": None,
            "pretrain_time_s": None,
            "finetune_time_s": None,
            "total_time_s": round(total_time, 1),
            "status": f"failed: {e}",
        }

    # Save per-task CSV
    df = pd.DataFrame([row])
    df.to_csv(csv_path, index=False)
    print(f"\n  Saved: {csv_path}")


# ============================================================
# CLI entry point
# ============================================================
if __name__ == "__main__":
    # Get config index from SLURM or CLI argument
    if "SLURM_ARRAY_TASK_ID" in os.environ:
        config_idx = int(os.environ["SLURM_ARRAY_TASK_ID"])
    elif len(sys.argv) > 1:
        config_idx = int(sys.argv[1])
    else:
        print(f"Usage: python {sys.argv[0]} <config_idx>")
        print(f"  or set SLURM_ARRAY_TASK_ID environment variable")
        print(f"\nTotal configs: {len(ALL_CONFIGS)}")
        print(f"\nRandom search space (seed={SEED}):")
        for k, v in SEARCH_SPACE.items():
            if isinstance(v, tuple):
                print(f"  {k}: uniform_int({v[0]}, {v[1]})")
            else:
                print(f"  {k}: choice({v})")
        print(f"\nConfig list:")
        for cfg in ALL_CONFIGS:
            baseline = " [A6 BASELINE]" if cfg["config_idx"] == 0 else ""
            print(f"  [{cfg['config_idx']:3d}] {cfg['name']}{baseline}")
        sys.exit(1)

    if config_idx < 0 or config_idx >= len(ALL_CONFIGS):
        print(f"ERROR: config_idx={config_idx} out of range [0, {len(ALL_CONFIGS)})")
        sys.exit(1)

    print(f"SLURM job: {os.environ.get('SLURM_JOB_ID', 'local')}")
    print(f"SLURM array task: {os.environ.get('SLURM_ARRAY_TASK_ID', 'N/A')}")
    print(f"GPU: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU'}")
    print(f"Config index: {config_idx}/{len(ALL_CONFIGS)}")
    print()

    run_one_config(config_idx)
