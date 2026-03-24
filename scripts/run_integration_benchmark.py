#!/usr/bin/env python
"""
run_integration_benchmark.py -- Generate integration benchmark reports from
pre-computed per-method h5ad files.

Reads ONLY obsm arrays via h5py (never loads the full X matrix) from
per-method h5ad files, assembles a lightweight AnnData with all embeddings,
subsamples, computes scib-style metrics, and generates an HTML report.

File naming convention: {sample}_{method}.h5ad
  e.g. robin_cell_harmony.h5ad, robin_cell_scvi.h5ad

Usage
-----
python scripts/run_integration_benchmark.py \
    --results-dir /path/to/preprocessing/ \
    --samples robin_008um robin_cell all_cell \
    --batch-key sample \
    --subsample-n 50000
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata import AnnData

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

# Arrow string coercion (needed when reading h5ad written on newer pandas)
pd.set_option("future.infer_string", False)
try:
    pd.options.mode.string_storage = "python"
except Exception:
    pass

# Method name -> (obsm embedding key, display name)
METHOD_CONFIG = {
    "harmony": ("X_pca_harmony", "Harmony"),
    "scvi": ("X_scVI", "scVI"),
    "resolvi": ("X_resolVI", "resolVI"),
}

# Unintegrated baseline key (present in all method files)
BASELINE_KEY = "X_pca"
BASELINE_NAME = "Unintegrated (PCA)"


def _read_obs_index(h5path: Path) -> np.ndarray:
    """Read the obs index (cell barcodes) from an h5ad file via h5py."""
    with h5py.File(h5path, "r") as f:
        obs_group = f["obs"]
        # _index stores the obs_names in h5ad format
        index_key = obs_group.attrs.get("_index", "_index")
        if isinstance(index_key, bytes):
            index_key = index_key.decode()
        if index_key in obs_group:
            idx = obs_group[index_key][:]
        elif "_index" in obs_group:
            idx = obs_group["_index"][:]
        else:
            raise KeyError(f"Cannot find obs index in {h5path}")
        # Decode bytes to str if needed
        if idx.dtype.kind == "S" or idx.dtype.kind == "O":
            idx = np.array([x.decode() if isinstance(x, bytes) else x for x in idx])
        return idx


def _read_obsm_array(h5path: Path, key: str) -> np.ndarray | None:
    """Read a single obsm array from an h5ad file via h5py. Returns None if missing."""
    with h5py.File(h5path, "r") as f:
        if "obsm" not in f:
            return None
        obsm = f["obsm"]
        if key not in obsm:
            return None
        return obsm[key][:]


def _read_obs_column(h5path: Path, col: str) -> np.ndarray | None:
    """Read a single obs column from an h5ad file via h5py. Returns None if missing."""
    with h5py.File(h5path, "r") as f:
        obs = f["obs"]
        if col not in obs:
            return None
        ds = obs[col]

        # Handle categorical encoding (h5ad stores categories + codes)
        if isinstance(ds, h5py.Group):
            # Categorical: has 'codes' and 'categories' datasets
            codes = ds["codes"][:]
            categories = ds["categories"][:]
            if categories.dtype.kind in ("S", "O"):
                categories = np.array(
                    [x.decode() if isinstance(x, bytes) else x for x in categories]
                )
            result = categories[codes]
            # Mask -1 codes (missing) as "NA"
            mask = codes < 0
            if mask.any():
                result = result.astype(object)
                result[mask] = "NA"
            return result
        else:
            data = ds[:]
            if data.dtype.kind in ("S", "O"):
                data = np.array([x.decode() if isinstance(x, bytes) else x for x in data])
            return data


def _find_base_adata(results_dir: Path, sample: str) -> Path | None:
    """Find the filtered base adata for a sample.

    Tries several naming conventions:
    - adata.{sample}.filtered.h5ad (sc_tools convention)
    - {sample}_harmony.h5ad (any method file works as a base)
    """
    # Look in parent directory and results directory
    for search_dir in [results_dir.parent, results_dir]:
        for pattern in [
            f"adata.{sample}.filtered.h5ad",
            f"{sample}.filtered.h5ad",
            f"adata_{sample}_filtered.h5ad",
        ]:
            candidate = search_dir / pattern
            if candidate.exists():
                return candidate

    # Fall back to any method file (we only read obs from it, not X)
    for method in METHOD_CONFIG:
        candidate = results_dir / f"{sample}_{method}.h5ad"
        if candidate.exists():
            return candidate

    return None


def _coerce_arrow_strings_df(df: pd.DataFrame) -> pd.DataFrame:
    """Convert Arrow-backed StringDtype columns to plain object dtype."""
    for col in df.columns:
        if hasattr(df[col].dtype, "name") and "string" in str(df[col].dtype).lower():
            df[col] = df[col].astype(object)
    return df


def load_embeddings_h5py(
    results_dir: Path,
    sample: str,
    batch_key: str,
) -> tuple[AnnData, dict[str, str]]:
    """Load obs metadata and obsm embeddings via h5py into a lightweight AnnData.

    Never loads the X matrix. Reads embeddings from per-method h5ad files
    and assembles them into a single AnnData.

    Returns
    -------
    tuple[AnnData, dict[str, str]]
        Lightweight AnnData (X is zeros) with all embeddings in obsm,
        and a dict mapping display name -> obsm key.
    """
    import anndata as ad

    # Find a base file to get obs metadata
    base_path = _find_base_adata(results_dir, sample)
    if base_path is None:
        raise FileNotFoundError(
            f"No base adata found for sample '{sample}' in {results_dir}. "
            f"Looked for adata.{sample}.filtered.h5ad and method files."
        )

    log.info("Using base file for obs metadata: %s", base_path)

    # Read obs index from base file
    obs_index = _read_obs_index(base_path)
    n_obs = len(obs_index)
    log.info("Base file has %d cells", n_obs)

    # Read batch key column
    batch_values = _read_obs_column(base_path, batch_key)
    if batch_values is None:
        raise KeyError(
            f"Batch key '{batch_key}' not found in obs of {base_path}. Check --batch-key argument."
        )

    # Build minimal obs DataFrame
    obs_df = pd.DataFrame({batch_key: batch_values}, index=obs_index)
    obs_df = _coerce_arrow_strings_df(obs_df)

    # Try to read additional useful columns (leiden, celltype, bio_label)
    # First try the base file, then fall back to method files for missing columns
    extra_cols_wanted = ["leiden", "celltype", "cell_type", "cluster", "bio_label"]
    for extra_col in extra_cols_wanted:
        vals = _read_obs_column(base_path, extra_col)
        if vals is not None:
            obs_df[extra_col] = vals

    # For any still-missing columns, try reading from method files
    missing_extras = [c for c in extra_cols_wanted if c not in obs_df.columns]
    if missing_extras:
        for method in METHOD_CONFIG:
            method_path = results_dir / f"{sample}_{method}.h5ad"
            if not method_path.exists():
                continue
            for col in list(missing_extras):
                vals = _read_obs_column(method_path, col)
                if vals is not None and len(vals) == n_obs:
                    obs_df[col] = vals
                    missing_extras.remove(col)
            if not missing_extras:
                break

    # Create lightweight AnnData (no X)
    adata = ad.AnnData(
        X=np.zeros((n_obs, 1), dtype=np.float32),
        obs=obs_df,
    )

    # Collect embeddings from method files
    embedding_keys: dict[str, str] = {}
    baseline_loaded = False

    for method, (emb_key, display_name) in METHOD_CONFIG.items():
        method_path = results_dir / f"{sample}_{method}.h5ad"
        if not method_path.exists():
            log.warning("Method file not found, skipping: %s", method_path)
            continue

        # Verify cell alignment
        method_index = _read_obs_index(method_path)
        if len(method_index) != n_obs:
            log.warning(
                "Cell count mismatch for %s: base=%d, method=%d. Skipping.",
                method,
                n_obs,
                len(method_index),
            )
            continue

        if not np.array_equal(method_index, obs_index):
            log.warning(
                "Cell index mismatch for %s (first 3: base=%s, method=%s). "
                "Proceeding assuming same order.",
                method,
                obs_index[:3],
                method_index[:3],
            )

        # Read the method embedding
        emb_array = _read_obsm_array(method_path, emb_key)
        if emb_array is not None:
            adata.obsm[emb_key] = emb_array
            embedding_keys[display_name] = emb_key
            log.info(
                "Loaded %s embedding: %s, shape=%s",
                display_name,
                emb_key,
                emb_array.shape,
            )
        else:
            log.warning("Embedding key '%s' not found in %s", emb_key, method_path)

        # Also grab the baseline PCA from the first method file if not yet loaded
        if not baseline_loaded:
            pca_array = _read_obsm_array(method_path, BASELINE_KEY)
            if pca_array is not None:
                adata.obsm[BASELINE_KEY] = pca_array
                embedding_keys[BASELINE_NAME] = BASELINE_KEY
                baseline_loaded = True
                log.info(
                    "Loaded baseline PCA from %s, shape=%s",
                    method_path.name,
                    pca_array.shape,
                )

    if not embedding_keys:
        raise RuntimeError(f"No embeddings loaded for sample '{sample}'")

    return adata, embedding_keys


def filter_nan_rows(
    embedding: np.ndarray,
    name: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Filter out rows with any NaN values from an embedding.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        (valid_mask, filtered_embedding)
    """
    nan_mask = np.isnan(embedding).any(axis=1)
    n_nan = nan_mask.sum()
    if n_nan > 0:
        log.info(
            "Embedding '%s': dropping %d / %d NaN rows (%.1f%%)",
            name,
            n_nan,
            len(embedding),
            100 * n_nan / len(embedding),
        )
    valid_mask = ~nan_mask
    return valid_mask, embedding[valid_mask]


def subsample_adata(
    adata: AnnData,
    batch_key: str,
    n: int = 50_000,
    seed: int = 42,
) -> AnnData:
    """Stratified subsample preserving batch proportions."""
    if adata.n_obs <= n:
        log.info("Dataset has %d cells (<= %d), no subsampling needed", adata.n_obs, n)
        return adata.copy()

    from sc_tools.bm.integration import _stratified_subsample

    return _stratified_subsample(adata, key=batch_key, n=n, seed=seed)


def run_benchmark_for_sample(
    results_dir: Path,
    sample: str,
    batch_key: str,
    figures_dir: Path,
    subsample_n: int,
    modality: str,
    bio_label: str | None = None,
) -> Path:
    """Run the integration benchmark for a single sample.

    Returns the path to the generated HTML report.
    """
    from sc_tools.qc.report import generate_post_integration_report

    log.info("=" * 60)
    log.info("Processing sample: %s", sample)
    log.info("=" * 60)

    # Step 1: Load embeddings via h5py
    adata, embedding_keys = load_embeddings_h5py(results_dir, sample, batch_key)
    log.info(
        "Loaded %d embeddings for %d cells: %s",
        len(embedding_keys),
        adata.n_obs,
        list(embedding_keys.keys()),
    )

    # Step 2: Subsample ONCE before metrics
    adata_sub = subsample_adata(adata, batch_key=batch_key, n=subsample_n)

    # Auto-detect bio label column
    if bio_label is None:
        for cand in ["bio_label", "celltype", "cell_type", "cluster"]:
            if cand in adata_sub.obs.columns:
                bio_label = cand
                log.info("Auto-detected bio label column: %s", bio_label)
                break
    if bio_label and bio_label not in adata_sub.obs.columns:
        log.warning("Bio label column '%s' not in obs, skipping bio metrics", bio_label)
        bio_label = None

    # Step 3: Filter NaN rows per embedding and compute metrics
    # We need to handle NaN rows (resolVI) per-embedding.
    # compare_integrations doesn't handle NaNs internally, so we compute
    # metrics manually for each embedding, filtering NaNs per-embedding.
    from sc_tools.bm.integration import compute_composite_score, compute_integration_metrics

    rows = []
    for display_name, emb_key in embedding_keys.items():
        if emb_key not in adata_sub.obsm:
            log.warning("Embedding %s not in subsampled adata, skipping", emb_key)
            continue

        emb = np.asarray(adata_sub.obsm[emb_key])
        valid_mask, emb_clean = filter_nan_rows(emb, display_name)

        if emb_clean.shape[0] < 100:
            log.warning(
                "Embedding '%s' has only %d valid cells after NaN filter, skipping",
                display_name,
                emb_clean.shape[0],
            )
            continue

        # Create a temporary adata slice with only valid cells for this embedding
        adata_tmp = adata_sub[valid_mask].copy()
        adata_tmp.obsm[emb_key] = emb_clean

        try:
            metrics = compute_integration_metrics(
                adata_tmp,
                emb_key,
                batch_key,
                celltype_key=bio_label,
                use_scib="auto",
            )
            composite = compute_composite_score(metrics)
            row = {"method": display_name, "embedding_key": emb_key}
            row.update(metrics)
            row.update(composite)
            row["n_cells_used"] = emb_clean.shape[0]
            rows.append(row)
            log.info(
                "  %s: overall=%.3f, batch=%.3f (n=%d cells)",
                display_name,
                composite["overall_score"],
                composite["batch_score"],
                emb_clean.shape[0],
            )
        except Exception:
            log.warning("Metrics computation failed for %s", display_name, exc_info=True)

    if not rows:
        raise RuntimeError(f"No metrics computed for sample '{sample}'")

    comparison_df = pd.DataFrame(rows)
    comparison_df = comparison_df.sort_values("overall_score", ascending=False).reset_index(
        drop=True
    )

    log.info(
        "Benchmark results:\n%s",
        comparison_df[["method", "overall_score", "batch_score"]].to_string(),
    )

    # Step 4: Generate report
    # Use the subsampled adata (with all embeddings) for the report
    sample_figures_dir = figures_dir / sample
    sample_figures_dir.mkdir(parents=True, exist_ok=True)

    report_path = generate_post_integration_report(
        adata_sub,
        output_dir=sample_figures_dir,
        embedding_keys=embedding_keys,
        batch_key=batch_key,
        celltype_key=bio_label,
        sample_col=batch_key,
        modality=modality,
        title=f"Integration Benchmark: {sample}",
        comparison_df=comparison_df,
    )

    log.info("Report saved: %s", report_path)

    # Also save the comparison CSV
    csv_path = sample_figures_dir / f"integration_benchmark_{sample}.csv"
    comparison_df.to_csv(csv_path, index=False)
    log.info("Metrics CSV saved: %s", csv_path)

    return report_path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument(
        "--results-dir",
        required=True,
        type=Path,
        help="Directory containing per-method h5ad files ({sample}_{method}.h5ad)",
    )
    p.add_argument(
        "--samples",
        required=True,
        nargs="+",
        help="Sample prefixes to process (e.g. robin_008um robin_cell all_cell)",
    )
    p.add_argument(
        "--batch-key",
        default="sample",
        help="obs column for batch correction (default: sample)",
    )
    p.add_argument(
        "--figures-dir",
        type=Path,
        default=None,
        help="Output directory for figures/reports (default: results-dir/../figures/)",
    )
    p.add_argument(
        "--subsample-n",
        type=int,
        default=50_000,
        help="Number of cells to subsample for metrics (default: 50000)",
    )
    p.add_argument(
        "--modality",
        default="visium_hd",
        help="Data modality (default: visium_hd)",
    )
    p.add_argument(
        "--bio-label",
        default=None,
        help="obs column for bio conservation metrics (default: auto-detect from bio_label, celltype, cell_type, cluster)",
    )
    return p.parse_args()


def main():
    args = parse_args()

    results_dir = args.results_dir.resolve()
    if not results_dir.exists():
        log.error("Results directory does not exist: %s", results_dir)
        sys.exit(1)

    # Default figures dir
    if args.figures_dir is not None:
        figures_dir = args.figures_dir.resolve()
    else:
        figures_dir = results_dir.parent / "figures"

    log.info("Results dir: %s", results_dir)
    log.info("Figures dir: %s", figures_dir)
    log.info("Samples: %s", args.samples)
    log.info("Batch key: %s", args.batch_key)
    log.info("Subsample N: %d", args.subsample_n)

    # Process each sample independently
    reports: list[Path] = []
    for sample in args.samples:
        try:
            report_path = run_benchmark_for_sample(
                results_dir=results_dir,
                sample=sample,
                batch_key=args.batch_key,
                figures_dir=figures_dir,
                subsample_n=args.subsample_n,
                modality=args.modality,
                bio_label=args.bio_label,
            )
            reports.append(report_path)
        except Exception:
            log.error("Failed to process sample '%s'", sample, exc_info=True)

    log.info("=" * 60)
    log.info("Completed %d / %d samples", len(reports), len(args.samples))
    for rp in reports:
        log.info("  Report: %s", rp)


if __name__ == "__main__":
    main()
