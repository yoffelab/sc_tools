#!/usr/bin/env python3
"""Aggregate scANVI hyperparameter search results.

Reads all config_*.csv files from results/scanvi_search/,
combines into a single ranked table, and identifies the optimal config.

Usage:
    python scripts/aggregate_scanvi_search.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

WORKDIR = Path("/home/fs01/juk4007/elementolab/projects/ibd_spatial")
SEARCH_DIR = WORKDIR / "results" / "scanvi_search"


def main():
    csv_files = sorted(SEARCH_DIR.glob("config_*.csv"))

    if not csv_files:
        print(f"ERROR: No config_*.csv files found in {SEARCH_DIR}")
        sys.exit(1)

    print(f"Found {len(csv_files)} result files")

    # Read and combine
    dfs = []
    for f in csv_files:
        try:
            df = pd.read_csv(f)
            dfs.append(df)
        except Exception as e:
            print(f"  WARNING: Could not read {f.name}: {e}")

    if not dfs:
        print("ERROR: No valid CSV files")
        sys.exit(1)

    results = pd.concat(dfs, ignore_index=True)
    print(f"Total configs: {len(results)}")

    # Split success/fail
    success = results[results["status"] == "success"].copy()
    failed = results[results["status"] != "success"].copy()
    print(f"Successful: {len(success)}, Failed: {len(failed)}")

    if len(failed) > 0:
        print("\nFailed configs:")
        for _, row in failed.iterrows():
            print(f"  [{row['config_idx']:2d}] {row['name']}: {row['status']}")

    if len(success) == 0:
        print("ERROR: No successful configs")
        sys.exit(1)

    # ============================================================
    # Rank by ct_broad_asw (primary objective)
    # ============================================================
    success = success.sort_values("ct_broad_asw", ascending=False).reset_index(drop=True)
    success["rank"] = range(1, len(success) + 1)

    print("\n" + "=" * 90)
    print("FULL RESULTS (ranked by ct_broad_asw)")
    print("=" * 90)

    display_cols = [
        "rank", "name", "n_layers", "n_latent", "n_hidden",
        "classification_ratio", "pretrain_epochs_max", "finetune_epochs_max",
        "ct_broad_asw", "ct_asw", "batch_score", "platform_entropy",
        "pred_accuracy", "total_time_s",
    ]
    print(success[display_cols].to_string(index=False))

    # ============================================================
    # Best config
    # ============================================================
    best = success.iloc[0]
    print(f"\n{'=' * 90}")
    print(f"OPTIMAL CONFIG: {best['name']}")
    print(f"{'=' * 90}")
    print(f"  n_layers             = {int(best['n_layers'])}")
    print(f"  n_latent             = {int(best['n_latent'])}")
    print(f"  n_hidden             = {int(best['n_hidden'])}")
    print(f"  classification_ratio = {int(best['classification_ratio'])}")
    print(f"  pretrain_epochs      = {int(best['pretrain_epochs_max'])}")
    print(f"  finetune_epochs      = {int(best['finetune_epochs_max'])}")
    print(f"  ---")
    print(f"  ct_broad_asw         = {best['ct_broad_asw']:.4f}")
    print(f"  ct_asw               = {best['ct_asw']:.4f}" if pd.notna(best.get("ct_asw")) else "  ct_asw               = N/A")
    print(f"  batch_score          = {best['batch_score']:.4f}")
    print(f"  platform_entropy     = {best['platform_entropy']:.4f}")
    print(f"  pred_accuracy        = {best['pred_accuracy']:.4f}" if pd.notna(best.get("pred_accuracy")) else "  pred_accuracy        = N/A")
    print(f"  total_time           = {best['total_time_s']:.0f}s")

    # Comparison to A6 baseline
    print(f"\n  vs A6 baseline (ct_broad_asw=0.189):")
    delta = best["ct_broad_asw"] - 0.189
    print(f"  Delta ct_broad_asw   = {delta:+.4f}")

    # ============================================================
    # Parameter sensitivity analysis
    # ============================================================
    print(f"\n{'=' * 90}")
    print("PARAMETER SENSITIVITY (mean ct_broad_asw by parameter value)")
    print("=" * 90)

    for param in ["n_layers", "n_latent", "n_hidden", "classification_ratio",
                  "pretrain_epochs_max", "finetune_epochs_max"]:
        grouped = success.groupby(param)["ct_broad_asw"].agg(["mean", "std", "count", "max"])
        print(f"\n  {param}:")
        for val, row in grouped.iterrows():
            bar = "#" * int(row["mean"] * 100)
            print(f"    {val:>4}: mean={row['mean']:.4f} +/- {row['std']:.4f}  "
                  f"max={row['max']:.4f}  n={int(row['count'])}  {bar}")

    # ============================================================
    # Top 5 vs Bottom 5
    # ============================================================
    print(f"\n{'=' * 90}")
    print("TOP 5 vs BOTTOM 5")
    print("=" * 90)

    top5 = success.head(5)
    bot5 = success.tail(5)

    print("\nTop 5:")
    for _, row in top5.iterrows():
        print(f"  [{int(row['config_idx']):2d}] {row['name']}: "
              f"ct_broad={row['ct_broad_asw']:.4f}, batch={row['batch_score']:.4f}, "
              f"entropy={row['platform_entropy']:.4f}")

    print("\nBottom 5:")
    for _, row in bot5.iterrows():
        print(f"  [{int(row['config_idx']):2d}] {row['name']}: "
              f"ct_broad={row['ct_broad_asw']:.4f}, batch={row['batch_score']:.4f}, "
              f"entropy={row['platform_entropy']:.4f}")

    # ============================================================
    # Pairwise correlations
    # ============================================================
    print(f"\n{'=' * 90}")
    print("METRIC CORRELATIONS")
    print("=" * 90)

    metric_cols = ["ct_broad_asw", "batch_score", "platform_entropy", "pred_accuracy"]
    available = [c for c in metric_cols if c in success.columns and success[c].notna().sum() > 5]
    if len(available) > 1:
        corr = success[available].corr()
        print(corr.round(3).to_string())

    # ============================================================
    # Save
    # ============================================================
    out_path = SEARCH_DIR / "scanvi_search_results.csv"
    success.to_csv(out_path, index=False)
    print(f"\nSaved combined results: {out_path}")

    # Also save a summary
    summary_path = SEARCH_DIR / "scanvi_search_summary.txt"
    with open(summary_path, "w") as f:
        f.write(f"scANVI Hyperparameter Search Summary\n")
        f.write(f"{'=' * 50}\n\n")
        f.write(f"Total configs: {len(results)}\n")
        f.write(f"Successful: {len(success)}\n")
        f.write(f"Failed: {len(failed)}\n\n")
        f.write(f"Optimal config: {best['name']}\n")
        f.write(f"  n_layers={int(best['n_layers'])}, n_latent={int(best['n_latent'])}, "
                f"n_hidden={int(best['n_hidden'])}\n")
        f.write(f"  classification_ratio={int(best['classification_ratio'])}\n")
        f.write(f"  pretrain={int(best['pretrain_epochs_max'])}, "
                f"finetune={int(best['finetune_epochs_max'])}\n\n")
        f.write(f"  ct_broad_asw = {best['ct_broad_asw']:.4f}\n")
        f.write(f"  batch_score  = {best['batch_score']:.4f}\n")
        f.write(f"  entropy      = {best['platform_entropy']:.4f}\n")
        f.write(f"  vs A6 delta  = {delta:+.4f}\n")
    print(f"Saved summary: {summary_path}")

    # Missing configs (if any tasks did not produce output)
    expected = set(range(64))
    found = set(success["config_idx"].astype(int).tolist()) | set(failed["config_idx"].astype(int).tolist())
    missing = expected - found
    if missing:
        print(f"\nWARNING: {len(missing)} configs have no output: {sorted(missing)}")
        print("These tasks may still be running or may have crashed without writing CSV.")


if __name__ == "__main__":
    main()
