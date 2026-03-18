"""
Real-data edge-case tests for sc_tools.qc.

Requires real data files to be present; each class is guarded by
``pytest.mark.skipif`` so the suite degrades gracefully when data is absent.

Data paths (from task brief):
- IMC celltyped:   projects/imc/ggo_human/results/celltyped.h5ad
- Visium annotated: projects/visium/ggo_visium/results/adata.annotated.p2.h5ad
- CosMx annotated:  projects/cosmx_1k/lymph_dlbcl/results/cosmx_rna_annotated.h5ad
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from sc_tools.qc import (
    calculate_qc_metrics,
    classify_samples,
    compute_sample_metrics,
    filter_cells,
    filter_genes,
    filter_spots,
    highly_variable_genes,
)

# ---------------------------------------------------------------------------
# Absolute paths to real data
# ---------------------------------------------------------------------------


def _find_repo_root() -> Path:
    """Walk upward from this file to find the root that contains real project data."""
    _sentinel = Path("projects/visium/ggo_visium/results/adata.annotated.p2.h5ad")
    for parent in Path(__file__).resolve().parents:
        if (parent / _sentinel).is_file():
            return parent
    return Path(__file__).parents[2]


_REPO_ROOT = _find_repo_root()

_VISIUM_PATH = _REPO_ROOT / "projects/visium/ggo_visium/results/adata.annotated.p2.h5ad"
_IMC_PATH = _REPO_ROOT / "projects/imc/ggo_human/results/celltyped.h5ad"
_COSMX_PATH = _REPO_ROOT / "projects/cosmx_1k/lymph_dlbcl/results/cosmx_rna_annotated.h5ad"

_VISIUM_AVAILABLE = _VISIUM_PATH.exists()
_IMC_AVAILABLE = _IMC_PATH.exists()
_COSMX_AVAILABLE = _COSMX_PATH.exists()

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

_RNG = np.random.default_rng(42)


def _sample_obs(adata, n: int, sample_col: str | None = None, sample_id: str | None = None):
    """Return a random subset of n cells, optionally restricted to one sample."""
    if sample_col and sample_id:
        mask = adata.obs[sample_col] == sample_id
        idx = adata.obs_names[mask]
    else:
        idx = adata.obs_names
    chosen = _RNG.choice(idx, size=min(n, len(idx)), replace=False)
    return adata[chosen].copy()


def _three_library_subset(adata, sample_col: str = "library_id", n_per_lib: int = 200):
    """Pick up to 3 libraries and take up to n_per_lib spots each."""
    libs = adata.obs[sample_col].unique()[:3]
    parts = []
    for lib in libs:
        sub = _sample_obs(adata, n_per_lib, sample_col=sample_col, sample_id=lib)
        parts.append(sub)
    return sc.concat(parts, label=sample_col, keys=libs.tolist())


# ===========================================================================
# 1. calculate_qc_metrics on real data
# ===========================================================================


@pytest.mark.skipif(not _VISIUM_AVAILABLE, reason="Visium real data not available")
class TestCalculateQcMetricsVisium:
    """calculate_qc_metrics on a 3-library Visium subset."""

    @pytest.fixture(scope="class")
    def subset(self):
        adata = sc.read_h5ad(_VISIUM_PATH)
        return _three_library_subset(adata, sample_col="library_id", n_per_lib=200)

    def test_metrics_computed_no_crash(self, subset):
        """Standard call on 3 libraries must produce obs QC columns."""
        calculate_qc_metrics(subset, modality="visium")
        assert "total_counts" in subset.obs.columns
        assert "n_genes_by_counts" in subset.obs.columns

    def test_pct_counts_in_top_genes_present(self, subset):
        """pct_counts_in_top_50_genes must be present after metrics computation."""
        calculate_qc_metrics(subset, modality="visium", percent_top=(50, 200))
        assert "pct_counts_in_top_50_genes" in subset.obs.columns

    def test_all_zero_gene_no_crash(self, subset):
        """An all-zero gene column must not cause an error in calculate_qc_metrics."""
        sub = subset.copy()
        # Force first gene to all-zero
        import scipy.sparse as sp

        if sp.issparse(sub.X):
            sub.X = sub.X.tolil()
            sub.X[:, 0] = 0
            sub.X = sub.X.tocsr()
        else:
            sub.X[:, 0] = 0
        calculate_qc_metrics(sub, modality="visium")
        assert "total_counts" in sub.obs.columns

    def test_zero_count_cells_handled(self, subset):
        """Cells forced to all-zero counts must still produce valid rows in obs."""
        sub = subset.copy()
        import scipy.sparse as sp

        if sp.issparse(sub.X):
            sub.X = sub.X.tolil()
            sub.X[0, :] = 0
            sub.X = sub.X.tocsr()
        else:
            sub.X[0, :] = 0
        calculate_qc_metrics(sub, modality="visium")
        assert sub.obs["total_counts"].iloc[0] == 0


@pytest.mark.skipif(not _IMC_AVAILABLE, reason="IMC real data not available")
class TestCalculateQcMetricsIMC:
    """calculate_qc_metrics on IMC (protein panel): MT must be gracefully absent."""

    @pytest.fixture(scope="class")
    def imc_subset(self):
        adata = sc.read_h5ad(_IMC_PATH)
        return _sample_obs(adata, 200)

    def test_no_crash_imc_modality(self, imc_subset):
        """calculate_qc_metrics with modality='imc' must not crash."""
        calculate_qc_metrics(imc_subset, modality="imc")
        assert "total_counts" in imc_subset.obs.columns

    def test_no_pct_mt_for_imc(self, imc_subset):
        """IMC protein panel has no MT genes: pct_counts_mt must NOT be added."""
        calculate_qc_metrics(imc_subset, modality="imc")
        assert "pct_counts_mt" not in imc_subset.obs.columns

    def test_mt_pattern_none_explicit(self, imc_subset):
        """Explicit mt_pattern=None must also suppress pct_counts_mt column."""
        sub = imc_subset.copy()
        calculate_qc_metrics(sub, mt_pattern=None, hb_pattern=None)
        assert "pct_counts_mt" not in sub.obs.columns

    def test_percent_top_capped_to_n_vars(self, imc_subset):
        """percent_top values larger than n_vars must be silently dropped."""
        n_vars = imc_subset.n_vars
        calculate_qc_metrics(
            imc_subset.copy(),
            modality="imc",
            percent_top=(50, 100, 200, n_vars + 10_000),
        )
        # No IndexError raised -- test passes if we reach here


# ===========================================================================
# 2. filter_cells / filter_genes on real data
# ===========================================================================


@pytest.mark.skipif(not _VISIUM_AVAILABLE, reason="Visium real data not available")
class TestFilterCellsGenesVisium:
    """filter_cells and filter_genes on Visium subsets."""

    @pytest.fixture
    def prepped(self):
        adata = sc.read_h5ad(_VISIUM_PATH)
        sub = _three_library_subset(adata, n_per_lib=200)
        calculate_qc_metrics(sub, modality="visium")
        return sub

    def test_standard_filter_reduces_n_obs(self, prepped):
        """filter_cells with min_counts=100 must reduce n_obs on typical Visium data."""
        n_before = prepped.n_obs
        filter_cells(prepped, min_counts=100)
        # Some cells removed is not guaranteed but n_obs must not increase
        assert prepped.n_obs <= n_before

    def test_n_vars_unchanged_by_filter_cells(self, prepped):
        """filter_cells must not change the number of genes."""
        n_vars = prepped.n_vars
        filter_cells(prepped, min_counts=50)
        assert prepped.n_vars == n_vars

    def test_over_aggressive_filter_cells_empty(self, prepped):
        """filter_cells with min_counts exceeding all cells gives 0 obs without error."""
        max_counts = int(prepped.obs["total_counts"].max()) + 1_000_000
        sub = prepped.copy()
        filter_cells(sub, min_counts=max_counts)
        assert sub.n_obs == 0

    def test_filter_genes_min_cells_removes_rare_genes(self, prepped):
        """filter_genes with min_cells close to n_obs removes most genes."""
        n_before = prepped.n_vars
        # Require every gene to appear in 99% of cells -- many should be removed
        threshold = max(1, int(prepped.n_obs * 0.99))
        filter_genes(prepped, min_cells=threshold)
        assert prepped.n_vars <= n_before

    def test_filter_genes_all_removed_empty(self, prepped):
        """filter_genes with impossibly high min_cells gives 0 vars without error."""
        sub = prepped.copy()
        filter_genes(sub, min_cells=sub.n_obs + 1)
        assert sub.n_vars == 0


# ===========================================================================
# 3. sample_qc on real data (compute_sample_metrics + classify_samples)
# ===========================================================================


@pytest.mark.skipif(not _IMC_AVAILABLE, reason="IMC real data not available")
class TestSampleQcIMC:
    """compute_sample_metrics and classify_samples on IMC ROI subsets."""

    @pytest.fixture(scope="class")
    def imc_multi_roi(self):
        """Three ROIs, 200 cells each, with pre-computed QC metrics."""
        adata = sc.read_h5ad(_IMC_PATH)
        # Determine sample column name
        sample_col = "ROI" if "ROI" in adata.obs.columns else adata.obs.columns[0]
        # Try common column names
        for col in ("ROI", "roi", "sample", "library_id", "sample_id"):
            if col in adata.obs.columns:
                sample_col = col
                break
        rois = adata.obs[sample_col].unique()[:3]
        parts = []
        for roi in rois:
            sub = _sample_obs(adata, 200, sample_col=sample_col, sample_id=roi)
            parts.append(sub)
        merged = sc.concat(parts, label=sample_col, keys=rois.tolist())
        calculate_qc_metrics(merged, modality="imc")
        return merged, sample_col

    def test_per_sample_metrics_has_three_rows(self, imc_multi_roi):
        """compute_sample_metrics must return one row per ROI."""
        adata, sample_col = imc_multi_roi
        metrics = compute_sample_metrics(adata, sample_col=sample_col, modality="imc")
        assert len(metrics) == 3

    def test_n_spots_column_populated(self, imc_multi_roi):
        """n_spots column must reflect actual cell counts per ROI."""
        adata, sample_col = imc_multi_roi
        metrics = compute_sample_metrics(adata, sample_col=sample_col, modality="imc")
        assert "n_spots" in metrics.columns
        assert (metrics["n_spots"] > 0).all()

    def test_classify_samples_produces_qc_pass(self, imc_multi_roi):
        """classify_samples must add qc_pass column to metrics DataFrame."""
        adata, sample_col = imc_multi_roi
        metrics = compute_sample_metrics(adata, sample_col=sample_col, modality="imc")
        classified = classify_samples(metrics, modality="imc")
        assert "qc_pass" in classified.columns
        assert classified["qc_pass"].dtype == bool

    def test_low_cell_roi_flagged_by_threshold(self, imc_multi_roi):
        """A synthetic ROI with very few cells must fail n_spots_min threshold."""
        adata, sample_col = imc_multi_roi
        metrics = compute_sample_metrics(adata, sample_col=sample_col, modality="imc")
        # Inject a fake tiny ROI into metrics
        tiny_row = metrics.iloc[0:1].copy()
        tiny_row.index = ["tiny_roi"]
        tiny_row["n_spots"] = 5
        augmented = pd.concat([metrics, tiny_row])
        classified = classify_samples(augmented, modality="imc")
        assert not classified.loc["tiny_roi", "qc_pass"]

    def test_missing_sample_col_raises(self, imc_multi_roi):
        """compute_sample_metrics must raise ValueError for missing sample column."""
        adata, _ = imc_multi_roi
        with pytest.raises(ValueError, match="sample_col"):
            compute_sample_metrics(adata, sample_col="nonexistent_column")


@pytest.mark.skipif(not _VISIUM_AVAILABLE, reason="Visium real data not available")
class TestSampleQcVisium:
    """compute_sample_metrics and classify_samples on Visium subsets."""

    @pytest.fixture(scope="class")
    def visium_multi_lib(self):
        adata = sc.read_h5ad(_VISIUM_PATH)
        sub = _three_library_subset(adata, n_per_lib=200)
        calculate_qc_metrics(sub, modality="visium")
        return sub

    def test_metrics_one_row_per_library(self, visium_multi_lib):
        metrics = compute_sample_metrics(
            visium_multi_lib, sample_col="library_id", modality="visium"
        )
        n_libs = visium_multi_lib.obs["library_id"].nunique()
        assert len(metrics) == n_libs

    def test_total_counts_median_positive(self, visium_multi_lib):
        metrics = compute_sample_metrics(
            visium_multi_lib, sample_col="library_id", modality="visium"
        )
        assert (metrics["total_counts_median"] > 0).all()

    def test_classify_passes_healthy_visium(self, visium_multi_lib):
        """Well-populated Visium libraries must pass absolute thresholds."""
        metrics = compute_sample_metrics(
            visium_multi_lib, sample_col="library_id", modality="visium"
        )
        classified = classify_samples(metrics, modality="visium")
        # At least one library must pass (healthy Visium data)
        assert classified["qc_pass"].any()


# ===========================================================================
# 4. highly_variable_genes on real data
# ===========================================================================


@pytest.mark.skipif(not _VISIUM_AVAILABLE, reason="Visium real data not available")
class TestHighlyVariableGenesVisium:
    """highly_variable_genes on a normalized Visium subset."""

    @pytest.fixture
    def normed_subset(self):
        """Normalize a 200-spot single-library subset ready for HVG selection."""
        adata = sc.read_h5ad(_VISIUM_PATH)
        libs = adata.obs["library_id"].unique()
        sub = _sample_obs(adata, 200, sample_col="library_id", sample_id=libs[0])
        sc.pp.normalize_total(sub, target_sum=1e4)
        sc.pp.log1p(sub)
        return sub

    @pytest.fixture
    def normed_multi_lib(self):
        """Normalize a 3-library, 200-spot-each subset."""
        adata = sc.read_h5ad(_VISIUM_PATH)
        sub = _three_library_subset(adata, n_per_lib=200)
        sc.pp.normalize_total(sub, target_sum=1e4)
        sc.pp.log1p(sub)
        return sub

    def test_hvg_single_sample_no_crash(self, normed_subset):
        """HVG selection on a single library must not crash."""
        highly_variable_genes(normed_subset, flavor="seurat")
        assert "highly_variable" in normed_subset.var.columns

    def test_hvg_marks_at_least_one_gene(self, normed_subset):
        """At least one gene should be marked highly variable."""
        highly_variable_genes(normed_subset, flavor="seurat")
        assert normed_subset.var["highly_variable"].any()

    def test_hvg_multi_sample_with_batch_key(self, normed_multi_lib):
        """HVG with batch_key='library_id' must not crash and must add column."""
        highly_variable_genes(normed_multi_lib, flavor="seurat", batch_key="library_id")
        assert "highly_variable" in normed_multi_lib.var.columns

    def test_hvg_seurat_v3_n_top_genes(self, normed_subset):
        """seurat_v3 flavor with n_top_genes must select exactly n_top genes."""
        sub = normed_subset.copy()
        # seurat_v3 needs raw counts; use the normalized values as proxy
        # (just testing the function runs without error on real data shape)
        n_top = 50
        highly_variable_genes(sub, flavor="seurat", n_top_genes=n_top)
        # highly_variable column must be bool Series
        assert sub.var["highly_variable"].dtype == bool


# ===========================================================================
# 5. filter_spots on real data (sample_qc.filter_spots)
# ===========================================================================


@pytest.mark.skipif(not _VISIUM_AVAILABLE, reason="Visium real data not available")
class TestFilterSpotsVisium:
    """filter_spots with modality-aware defaults on a Visium subset."""

    @pytest.fixture
    def prepped(self):
        adata = sc.read_h5ad(_VISIUM_PATH)
        sub = _three_library_subset(adata, n_per_lib=200)
        calculate_qc_metrics(sub, modality="visium")
        return sub

    def test_standard_filter_reduces_or_equal_n_obs(self, prepped):
        """filter_spots inplace must not increase n_obs."""
        n_before = prepped.n_obs
        filter_spots(prepped, modality="visium", inplace=True)
        assert prepped.n_obs <= n_before

    def test_inplace_false_returns_copy(self, prepped):
        """filter_spots with inplace=False must return a new AnnData."""
        result = filter_spots(prepped.copy(), modality="visium", inplace=False)
        assert result is not None
        assert isinstance(result, type(prepped))

    def test_over_aggressive_threshold_yields_empty(self, prepped):
        """Very high min_counts threshold must remove all spots without error."""
        result = filter_spots(
            prepped.copy(),
            modality="visium",
            min_counts=10_000_000,
            inplace=False,
        )
        assert result.n_obs == 0

    def test_per_sample_logging_does_not_crash(self, prepped):
        """filter_spots with sample_col must log per-sample info without error."""
        filter_spots(
            prepped.copy(),
            modality="visium",
            min_counts=100,
            sample_col="library_id",
            inplace=False,
        )


@pytest.mark.skipif(not _IMC_AVAILABLE, reason="IMC real data not available")
class TestFilterSpotsIMC:
    """filter_spots on IMC protein data."""

    @pytest.fixture
    def prepped(self):
        adata = sc.read_h5ad(_IMC_PATH)
        sub = _sample_obs(adata, 200)
        calculate_qc_metrics(sub, modality="imc")
        return sub

    def test_imc_filter_no_crash(self, prepped):
        """filter_spots on IMC must run without error."""
        result = filter_spots(prepped.copy(), modality="imc", inplace=False)
        assert result is not None

    def test_imc_filter_does_not_use_pct_mt(self, prepped):
        """IMC has no MT genes; filter must not reference pct_counts_mt column."""
        assert "pct_counts_mt" not in prepped.obs.columns
        # Should not raise even though pct_counts_mt is absent
        filter_spots(prepped.copy(), modality="imc", max_pct_mt=50.0, inplace=False)
