"""Phase 1 (qc_filter): QC filtering for robin cell segmentation data.

Reads results/adata.filtered.h5ad, applies xenium-like QC thresholds (single-cell),
generates pre-filter QC report, saves filtered results/adata.filtered.h5ad.

Run from project root or via SLURM.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
REPO_ROOT = PROJECT_ROOT.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def main():
    import scanpy as sc

    from sc_tools.qc import calculate_qc_metrics, filter_genes
    from sc_tools.qc.report import generate_pre_filter_report
    from sc_tools.qc.sample_qc import apply_qc_filter, classify_samples, compute_sample_metrics

    raw_path = PROJECT_ROOT / "results" / "adata.filtered.h5ad"
    if not raw_path.exists():
        logger.error("adata.filtered.h5ad not found at %s", raw_path)
        sys.exit(1)

    logger.info("Reading %s", raw_path)
    adata = sc.read_h5ad(raw_path)
    adata.obs_names_make_unique()
    logger.info("Loaded: %d cells x %d genes", adata.n_obs, adata.n_vars)

    # QC metrics
    calculate_qc_metrics(adata, mt_pattern="^MT-", hb_pattern="^HB[AB]")

    # Sample-level QC
    sample_col = "sample"
    metrics = compute_sample_metrics(adata, sample_col=sample_col, modality="xenium")
    classified = classify_samples(metrics, modality="xenium")
    n_fail = (classified["qc_pass"] == False).sum()  # noqa: E712
    logger.info("Sample QC: %d pass, %d fail", len(classified) - n_fail, n_fail)

    # Pre-filter QC report
    figures_dir = PROJECT_ROOT / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    try:
        report_path = generate_pre_filter_report(
            adata,
            metrics=metrics,
            classified=classified,
            output_dir=str(figures_dir / "QC"),
            sample_col=sample_col,
            modality="xenium",
        )
        logger.info("Pre-filter QC report: %s", report_path)
    except Exception as exc:
        logger.warning("Pre-filter report failed (non-fatal): %s", exc)

    # Apply QC filter (xenium-like thresholds for single-cell data)
    adata = apply_qc_filter(
        adata,
        classified=classified,
        sample_col=sample_col,
        modality="xenium",
        min_genes=5,
        min_counts=10,
    )
    logger.info("After QC filter: %d cells x %d genes", adata.n_obs, adata.n_vars)

    # Filter genes detected in very few cells
    filter_genes(adata, min_cells=10)
    logger.info("After gene filter: %d cells x %d genes", adata.n_obs, adata.n_vars)

    # Save (overwrite raw with filtered)
    adata.write_h5ad(raw_path)
    logger.info("Saved filtered adata: %s", raw_path)


if __name__ == "__main__":
    main()
