"""
Sample-level QC: spot filtering, per-sample metrics, pass/fail classification.

Provides:
- filter_spots: Remove low-quality spots/cells with modality-aware defaults.
- compute_sample_metrics: Per-sample aggregate QC metrics.
- classify_samples: Absolute threshold + MAD-based outlier detection.
- save_pass_fail_lists: Write pass/fail CSVs.
- apply_qc_filter: Full pipeline — backup, spot filter, sample removal, save.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData

__all__ = [
    "filter_spots",
    "compute_sample_metrics",
    "classify_samples",
    "save_pass_fail_lists",
    "apply_qc_filter",
]

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Default thresholds — intentionally very lenient
# ---------------------------------------------------------------------------

_SPOT_FILTER_DEFAULTS: dict[str, dict[str, Any]] = {
    "visium": {"min_counts": 50, "min_genes": 20, "max_pct_mt": None},
    "visium_hd": {"min_counts": 10, "min_genes": 5, "max_pct_mt": None},
    "visium_hd_cell": {"min_counts": 5, "min_genes": 3, "max_pct_mt": None},
    "xenium": {"min_counts": 5, "min_genes": 3, "max_pct_mt": None},
    "cosmx": {"min_counts": 5, "min_genes": 3, "max_pct_mt": None},
    "imc": {"min_counts": 1, "min_genes": 1, "max_pct_mt": None},
}

_SAMPLE_THRESHOLDS: dict[str, dict[str, Any]] = {
    "visium": {
        "n_genes_median_min": 50,
        "total_counts_median_min": 100,
        "pct_mt_median_max": 50.0,
        "n_spots_min": 50,
    },
    "visium_hd": {
        "n_genes_median_min": 5,
        "total_counts_median_min": 5,
        "pct_mt_median_max": 50.0,
        "n_spots_min": 500,
    },
    "visium_hd_cell": {
        "n_genes_median_min": 5,
        "total_counts_median_min": 10,
        "pct_mt_median_max": None,
        "n_spots_min": 50,
    },
    "xenium": {
        "n_genes_median_min": 5,
        "total_counts_median_min": 10,
        "pct_mt_median_max": None,
        "n_spots_min": 50,
    },
    "cosmx": {
        "n_genes_median_min": 5,
        "total_counts_median_min": 10,
        "pct_mt_median_max": None,
        "n_spots_min": 50,
    },
    "imc": {
        "n_genes_median_min": 3,
        "total_counts_median_min": 5,
        "pct_mt_median_max": None,
        "n_spots_min": 20,
    },
}


def _adaptive_mad_multiplier(n_samples: int, base_multiplier: float = 3.0) -> float:
    """Return conservative MAD multiplier scaled by cohort size."""
    if n_samples < 10:
        return 5.0
    elif n_samples < 20:
        return 4.0
    elif n_samples < 40:
        return base_multiplier
    else:
        return max(2.5, base_multiplier - 0.5)


def _mad(x: np.ndarray) -> float:
    """Median absolute deviation."""
    med = np.nanmedian(x)
    return float(np.nanmedian(np.abs(x - med)))


# ---------------------------------------------------------------------------
# filter_spots
# ---------------------------------------------------------------------------


def filter_spots(
    adata: AnnData,
    modality: str = "visium",
    min_counts: int | None = None,
    min_genes: int | None = None,
    max_pct_mt: float | None = None,
    sample_col: str | None = None,
    inplace: bool = True,
) -> AnnData | None:
    """
    Remove low-quality spots/cells using modality-aware defaults.

    Parameters
    ----------
    adata : AnnData
        Must have obs columns ``total_counts`` and ``n_genes_by_counts``
        (from ``calculate_qc_metrics``).
    modality : str
        One of ``visium``, ``visium_hd``, ``xenium``, ``cosmx``, ``imc``.
    min_counts, min_genes, max_pct_mt : int/float or None
        Override modality defaults. ``None`` means use the modality default
        (which itself may be ``None`` for max_pct_mt).
    sample_col : str or None
        If provided, log removal counts per sample.
    inplace : bool
        If True, filter in place and return None.

    Returns
    -------
    AnnData or None
        Filtered copy if ``inplace=False``; else None.
    """
    defaults = _SPOT_FILTER_DEFAULTS.get(modality, _SPOT_FILTER_DEFAULTS["visium"])
    mc = min_counts if min_counts is not None else defaults["min_counts"]
    mg = min_genes if min_genes is not None else defaults["min_genes"]
    mp = max_pct_mt if max_pct_mt is not None else defaults["max_pct_mt"]

    if not inplace:
        adata = adata.copy()

    n_before = adata.n_obs

    # Build boolean mask of spots to keep
    keep = np.ones(adata.n_obs, dtype=bool)

    if mc is not None and "total_counts" in adata.obs.columns:
        keep &= adata.obs["total_counts"].values >= mc

    if mg is not None and "n_genes_by_counts" in adata.obs.columns:
        keep &= adata.obs["n_genes_by_counts"].values >= mg

    if mp is not None and "pct_counts_mt" in adata.obs.columns:
        keep &= adata.obs["pct_counts_mt"].values <= mp

    from .report_utils import get_modality_terms

    terms = get_modality_terms(modality)
    _obs_label = terms["observations_lower"]

    if sample_col and sample_col in adata.obs.columns:
        removed = ~keep
        for sample in adata.obs[sample_col].unique():
            mask = adata.obs[sample_col] == sample
            n_rm = int(removed[mask].sum())
            if n_rm > 0:
                logger.info(
                    "filter_spots: %s — removed %d / %d %s",
                    sample,
                    n_rm,
                    int(mask.sum()),
                    _obs_label,
                )

    # Apply filter
    if not keep.all():
        adata._inplace_subset_obs(keep)

    n_after = adata.n_obs
    logger.info(
        "filter_spots (%s): %d -> %d %s (removed %d)",
        modality,
        n_before,
        n_after,
        _obs_label,
        n_before - n_after,
    )

    if not inplace:
        return adata
    return None


# ---------------------------------------------------------------------------
# compute_sample_metrics
# ---------------------------------------------------------------------------


def compute_sample_metrics(
    adata: AnnData,
    sample_col: str = "library_id",
    modality: str = "visium",
    spaceranger_dirs: dict[str, str | Path] | None = None,
) -> pd.DataFrame:
    """
    Compute per-sample aggregate QC metrics.

    Parameters
    ----------
    adata : AnnData
        Must have QC columns in obs (``total_counts``, ``n_genes_by_counts``,
        optionally ``pct_counts_mt``).
    sample_col : str
        Column in ``adata.obs`` identifying samples (default ``library_id``).
    modality : str
        Modality name (for future modality-specific metrics).
    spaceranger_dirs : dict or None
        Mapping sample name -> Space Ranger ``outs/`` directory. If provided,
        sequencing metrics are parsed from ``metrics_summary.csv``.

    Returns
    -------
    pd.DataFrame
        Indexed by sample with aggregate metric columns.
    """
    if sample_col not in adata.obs.columns:
        raise ValueError(f"sample_col={sample_col!r} not in adata.obs.columns")

    grouped = adata.obs.groupby(sample_col, observed=True)
    records = []

    for sample, grp in grouped:
        rec: dict[str, Any] = {"sample": sample, "n_spots": len(grp)}

        if "total_counts" in grp.columns:
            tc = grp["total_counts"].values.astype(float)
            rec["total_counts_median"] = float(np.nanmedian(tc))
            rec["total_counts_mean"] = float(np.nanmean(tc))
            rec["total_counts_sum"] = float(np.nansum(tc))

        if "n_genes_by_counts" in grp.columns:
            ng = grp["n_genes_by_counts"].values.astype(float)
            rec["n_genes_median"] = float(np.nanmedian(ng))
            rec["n_genes_mean"] = float(np.nanmean(ng))

        # Count unique detected genes for this sample
        sample_mask = adata.obs[sample_col] == sample
        sub_x = adata[sample_mask].X
        if hasattr(sub_x, "toarray"):
            sub_x = sub_x.toarray()
        rec["n_genes_detected"] = int(np.sum(np.asarray(sub_x).sum(axis=0) > 0))

        if "pct_counts_mt" in grp.columns:
            pmt = grp["pct_counts_mt"].values.astype(float)
            rec["pct_mt_median"] = float(np.nanmedian(pmt))
            rec["pct_mt_mean"] = float(np.nanmean(pmt))
            rec["pct_mt_max"] = float(np.nanmax(pmt))
            rec["pct_mt_gt5"] = float(np.nanmean(pmt > 5.0))
            rec["pct_mt_gt20"] = float(np.nanmean(pmt > 20.0))

        records.append(rec)

    metrics = pd.DataFrame(records).set_index("sample")
    metrics.index.name = sample_col

    # Parse Space Ranger metrics if provided
    if spaceranger_dirs:
        for sample, sr_dir in spaceranger_dirs.items():
            csv_path = Path(sr_dir) / "metrics_summary.csv"
            if sample in metrics.index and csv_path.exists():
                try:
                    sr = pd.read_csv(csv_path, thousands=",")
                    sr_dict = sr.iloc[0].to_dict()
                    for col in [
                        "Sequencing Saturation",
                        "Valid Barcodes",
                        "Reads Mapped Confidently to Genome",
                        "Mean Reads per Spot",
                        "Q30 Bases in Barcode",
                        "Q30 Bases in Probe Read",
                    ]:
                        if col in sr_dict:
                            val = sr_dict[col]
                            if isinstance(val, str):
                                val = val.replace("%", "").replace(",", "")
                            clean_col = col.lower().replace(" ", "_")
                            metrics.loc[sample, clean_col] = float(val)
                except Exception:
                    logger.warning("Could not parse Space Ranger metrics for %s", sample)

    return metrics


# ---------------------------------------------------------------------------
# classify_samples
# ---------------------------------------------------------------------------


def classify_samples(
    metrics: pd.DataFrame,
    modality: str = "visium",
    thresholds: dict[str, Any] | None = None,
    mad_multiplier: float = 3.0,
    min_cohort_size_for_outlier: int = 5,
) -> pd.DataFrame:
    """
    Classify samples as pass/fail using absolute thresholds and MAD outlier detection.

    Parameters
    ----------
    metrics : pd.DataFrame
        Output of ``compute_sample_metrics``.
    modality : str
        Modality for default thresholds.
    thresholds : dict or None
        Override default absolute thresholds. Keys match ``_SAMPLE_THRESHOLDS``.
    mad_multiplier : float
        Base MAD multiplier (adapted by cohort size).
    min_cohort_size_for_outlier : int
        Skip outlier detection if fewer samples than this.

    Returns
    -------
    pd.DataFrame
        Input DataFrame with added columns: ``qc_pass``, ``qc_fail_reasons``,
        ``qc_flag_absolute``, ``qc_flag_outlier``.
    """
    result = metrics.copy()
    n_samples = len(result)

    # Merge thresholds
    defaults = _SAMPLE_THRESHOLDS.get(modality, _SAMPLE_THRESHOLDS["visium"]).copy()
    if thresholds:
        defaults.update(thresholds)

    abs_flags = pd.Series(False, index=result.index)
    abs_reasons: dict[str, list[str]] = {s: [] for s in result.index}

    # Absolute thresholds
    if defaults.get("n_genes_median_min") is not None and "n_genes_median" in result.columns:
        bad = result["n_genes_median"] < defaults["n_genes_median_min"]
        for s in result.index[bad]:
            abs_reasons[s].append(
                f"n_genes_median={result.loc[s, 'n_genes_median']:.0f}"
                f" < {defaults['n_genes_median_min']}"
            )
        abs_flags |= bad

    if (
        defaults.get("total_counts_median_min") is not None
        and "total_counts_median" in result.columns
    ):
        bad = result["total_counts_median"] < defaults["total_counts_median_min"]
        for s in result.index[bad]:
            abs_reasons[s].append(
                f"total_counts_median={result.loc[s, 'total_counts_median']:.0f}"
                f" < {defaults['total_counts_median_min']}"
            )
        abs_flags |= bad

    if defaults.get("pct_mt_median_max") is not None and "pct_mt_median" in result.columns:
        bad = result["pct_mt_median"] > defaults["pct_mt_median_max"]
        for s in result.index[bad]:
            abs_reasons[s].append(
                f"pct_mt_median={result.loc[s, 'pct_mt_median']:.1f}%"
                f" > {defaults['pct_mt_median_max']}%"
            )
        abs_flags |= bad

    if defaults.get("n_spots_min") is not None and "n_spots" in result.columns:
        bad = result["n_spots"] < defaults["n_spots_min"]
        for s in result.index[bad]:
            abs_reasons[s].append(f"n_spots={result.loc[s, 'n_spots']} < {defaults['n_spots_min']}")
        abs_flags |= bad

    # Outlier detection (MAD-based)
    outlier_flags = pd.Series(False, index=result.index)
    outlier_reasons: dict[str, list[str]] = {s: [] for s in result.index}

    if n_samples >= min_cohort_size_for_outlier:
        effective_mad = _adaptive_mad_multiplier(n_samples, mad_multiplier)

        # Low outliers: n_genes_median, total_counts_median, n_spots
        for col in ["n_genes_median", "total_counts_median", "n_spots"]:
            if col not in result.columns:
                continue
            vals = result[col].values.astype(float)
            med = np.nanmedian(vals)
            mad_val = _mad(vals)
            if mad_val > 0:
                lower = med - effective_mad * mad_val
                bad = vals < lower
                for i, s in enumerate(result.index):
                    if bad[i]:
                        outlier_reasons[s].append(
                            f"{col}={vals[i]:.1f} < {lower:.1f} (MAD outlier)"
                        )
                outlier_flags |= pd.Series(bad, index=result.index)

        # High outlier: pct_mt_median
        if "pct_mt_median" in result.columns:
            vals = result["pct_mt_median"].values.astype(float)
            med = np.nanmedian(vals)
            mad_val = _mad(vals)
            if mad_val > 0:
                upper = med + effective_mad * mad_val
                bad = vals > upper
                for i, s in enumerate(result.index):
                    if bad[i]:
                        outlier_reasons[s].append(
                            f"pct_mt_median={vals[i]:.1f}% > {upper:.1f}% (MAD outlier)"
                        )
                outlier_flags |= pd.Series(bad, index=result.index)

    result["qc_flag_absolute"] = abs_flags.values
    result["qc_flag_outlier"] = outlier_flags.values
    result["qc_pass"] = ~(abs_flags | outlier_flags).values

    # Combine reasons
    all_reasons = []
    for s in result.index:
        reasons = abs_reasons[s] + outlier_reasons[s]
        all_reasons.append("; ".join(reasons))
    result["qc_fail_reasons"] = all_reasons

    return result


# ---------------------------------------------------------------------------
# save_pass_fail_lists
# ---------------------------------------------------------------------------


def save_pass_fail_lists(
    classified: pd.DataFrame,
    output_dir: str | Path,
    sample_col: str = "library_id",
) -> tuple[Path, Path]:
    """
    Write ``qc_sample_pass.csv`` and ``qc_sample_fail.csv``.

    Parameters
    ----------
    classified : pd.DataFrame
        Output of ``classify_samples``.
    output_dir : str or Path
        Directory for output CSVs.
    sample_col : str
        Name for the sample index column in output.

    Returns
    -------
    tuple of Path
        (pass_path, fail_path)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    passed = classified[classified["qc_pass"]].copy()
    failed = classified[~classified["qc_pass"]].copy()

    pass_path = output_dir / "qc_sample_pass.csv"
    fail_path = output_dir / "qc_sample_fail.csv"

    pass_cols = [c for c in passed.columns if c != "qc_fail_reasons"]
    passed[pass_cols].to_csv(pass_path)

    failed.to_csv(fail_path)

    logger.info(
        "QC classification: %d passed, %d failed -> %s",
        len(passed),
        len(failed),
        output_dir,
    )
    return pass_path, fail_path


# ---------------------------------------------------------------------------
# apply_qc_filter
# ---------------------------------------------------------------------------


def apply_qc_filter(
    adata: AnnData,
    classified: pd.DataFrame,
    sample_col: str = "library_id",
    modality: str = "visium",
    output_path: str | Path | None = None,
    backup_path: str | Path | None = None,
    min_counts: int | None = None,
    min_genes: int | None = None,
    max_pct_mt: float | None = None,
) -> AnnData:
    """
    Full QC pipeline: backup, spot-level filter, sample removal, save.

    Parameters
    ----------
    adata : AnnData
        Raw AnnData (will be modified in place).
    classified : pd.DataFrame
        Output of ``classify_samples`` with ``qc_pass`` column.
    sample_col : str
        Column in ``adata.obs`` identifying samples.
    modality : str
        Modality for spot-filter defaults.
    output_path : str or Path or None
        Save filtered AnnData here.
    backup_path : str or Path or None
        Save full (unfiltered) backup here before filtering.
    min_counts, min_genes, max_pct_mt
        Override spot-filter defaults.

    Returns
    -------
    AnnData
        Filtered AnnData (spots filtered, failed samples removed).
    """
    # 1. Save backup
    if backup_path is not None:
        bp = Path(backup_path)
        bp.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(bp)
        logger.info("Backup saved: %s (%d obs)", bp, adata.n_obs)

    # 2. Spot-level filtering
    filter_spots(
        adata,
        modality=modality,
        min_counts=min_counts,
        min_genes=min_genes,
        max_pct_mt=max_pct_mt,
        sample_col=sample_col,
        inplace=True,
    )

    # 3. Remove failed samples
    if sample_col in adata.obs.columns:
        failed = classified.index[~classified["qc_pass"]].tolist()
        if failed:
            n_before = adata.n_obs
            keep = ~adata.obs[sample_col].isin(failed)
            adata = adata[keep].copy()
            logger.info(
                "Removed %d failed samples (%d -> %d obs): %s",
                len(failed),
                n_before,
                adata.n_obs,
                failed,
            )

    # 4. Save filtered
    if output_path is not None:
        op = Path(output_path)
        op.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(op)
        logger.info("Filtered AnnData saved: %s (%d obs)", op, adata.n_obs)

    return adata
