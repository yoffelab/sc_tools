"""Tests for sc_tools.bm.integration — batch correction quality metrics."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from sc_tools.bm.integration import (
    compare_integrations,
    compute_composite_score,
    compute_integration_metrics,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _make_batched_adata(
    n_obs: int = 200,
    n_vars: int = 50,
    n_batches: int = 3,
    n_celltypes: int = 4,
) -> AnnData:
    """Create a synthetic AnnData with batch and celltype labels plus embeddings."""
    rng = np.random.RandomState(42)

    X = rng.randn(n_obs, n_vars).astype(np.float32)
    adata = AnnData(X=X)

    # Assign batches and cell types
    adata.obs["batch"] = [f"batch_{i % n_batches}" for i in range(n_obs)]
    adata.obs["celltype"] = [f"type_{i % n_celltypes}" for i in range(n_obs)]

    # "Well-mixed" embedding: cell types separate, batches mixed
    embedding_good = np.zeros((n_obs, 10), dtype=np.float32)
    for i in range(n_obs):
        ct = i % n_celltypes
        embedding_good[i] = rng.randn(10) * 0.5 + ct * 3  # cluster by celltype
    adata.obsm["X_good"] = embedding_good

    # "Poorly-mixed" embedding: batches separate (bad integration)
    embedding_bad = np.zeros((n_obs, 10), dtype=np.float32)
    for i in range(n_obs):
        bt = i % n_batches
        embedding_bad[i] = rng.randn(10) * 0.5 + bt * 5  # cluster by batch
    adata.obsm["X_bad"] = embedding_bad

    # PCA (unintegrated baseline)
    adata.obsm["X_pca"] = rng.randn(n_obs, 10).astype(np.float32)

    return adata


# ---------------------------------------------------------------------------
# compute_integration_metrics
# ---------------------------------------------------------------------------


class TestComputeIntegrationMetrics:
    def test_returns_dict(self):
        adata = _make_batched_adata()
        metrics = compute_integration_metrics(
            adata, "X_good", "batch", "celltype", use_scib="sklearn"
        )
        assert isinstance(metrics, dict)
        assert "asw_batch" in metrics
        assert "asw_celltype" in metrics
        assert "ari" in metrics
        assert "nmi" in metrics
        assert "pcr" in metrics
        assert "graph_connectivity" in metrics

    def test_metrics_in_0_1(self):
        adata = _make_batched_adata()
        metrics = compute_integration_metrics(
            adata, "X_good", "batch", "celltype", use_scib="sklearn"
        )
        for key, val in metrics.items():
            assert 0.0 <= val <= 1.0, f"{key}={val} not in [0,1]"

    def test_missing_embedding_raises(self):
        adata = _make_batched_adata()
        with pytest.raises(KeyError, match="X_nonexistent"):
            compute_integration_metrics(
                adata, "X_nonexistent", "batch", "celltype", use_scib="sklearn"
            )

    def test_missing_batch_key_raises(self):
        adata = _make_batched_adata()
        with pytest.raises(KeyError, match="missing_batch"):
            compute_integration_metrics(
                adata, "X_good", "missing_batch", "celltype", use_scib="sklearn"
            )

    def test_missing_celltype_key_skips_bio(self):
        """When celltype_key is not in obs, bio metrics are skipped (not raised)."""
        adata = _make_batched_adata()
        metrics = compute_integration_metrics(
            adata, "X_good", "batch", "missing_ct", use_scib="sklearn"
        )
        # Should return batch-only metrics (missing_ct not in obs -> treated as absent)
        assert "asw_batch" in metrics
        assert "pcr" in metrics
        assert "asw_celltype" not in metrics

    def test_well_mixed_better_asw_batch(self):
        adata = _make_batched_adata(n_obs=300)
        good = compute_integration_metrics(adata, "X_good", "batch", "celltype", use_scib="sklearn")
        bad = compute_integration_metrics(adata, "X_bad", "batch", "celltype", use_scib="sklearn")
        # Well-mixed should have better (higher) ASW batch score
        assert good["asw_batch"] >= bad["asw_batch"] - 0.1, (
            f"Expected good ASW batch >= bad, got {good['asw_batch']} vs {bad['asw_batch']}"
        )

    def test_well_mixed_better_bio_ari(self):
        adata = _make_batched_adata(n_obs=300)
        good = compute_integration_metrics(adata, "X_good", "batch", "celltype", use_scib="sklearn")
        bad = compute_integration_metrics(adata, "X_bad", "batch", "celltype", use_scib="sklearn")
        # Well-mixed (celltype-separated) should have better ARI
        assert good["ari"] >= bad["ari"] - 0.1, (
            f"Expected good ARI >= bad, got {good['ari']} vs {bad['ari']}"
        )


# ---------------------------------------------------------------------------
# compute_composite_score
# ---------------------------------------------------------------------------


class TestCompositeScore:
    def test_returns_keys(self):
        metrics = {
            "asw_batch": 0.8,
            "pcr": 0.7,
            "graph_connectivity": 0.9,
            "asw_celltype": 0.85,
            "ari": 0.7,
            "nmi": 0.75,
        }
        result = compute_composite_score(metrics)
        assert "batch_score" in result
        assert "bio_score" in result
        assert "overall_score" in result

    def test_overall_in_0_1(self):
        metrics = {
            "asw_batch": 0.5,
            "pcr": 0.5,
            "asw_celltype": 0.5,
            "ari": 0.5,
            "nmi": 0.5,
        }
        result = compute_composite_score(metrics)
        assert 0.0 <= result["overall_score"] <= 1.0

    def test_weights_affect_score(self):
        metrics = {
            "asw_batch": 1.0,
            "pcr": 1.0,
            "asw_celltype": 0.0,
            "ari": 0.0,
            "nmi": 0.0,
        }
        # High batch weight should favor batch metrics
        high_batch = compute_composite_score(metrics, batch_weight=0.9, bio_weight=0.1)
        high_bio = compute_composite_score(metrics, batch_weight=0.1, bio_weight=0.9)
        assert high_batch["overall_score"] > high_bio["overall_score"]

    def test_empty_metrics(self):
        result = compute_composite_score({})
        assert result["overall_score"] == 0.0

    def test_partial_metrics(self):
        # Only batch metrics
        result = compute_composite_score({"asw_batch": 0.8, "pcr": 0.6})
        assert result["batch_score"] > 0
        assert result["bio_score"] == 0.0


# ---------------------------------------------------------------------------
# compare_integrations
# ---------------------------------------------------------------------------


class TestCompareIntegrations:
    def test_returns_sorted_df(self):
        adata = _make_batched_adata()
        df = compare_integrations(
            adata,
            {"good": "X_good", "bad": "X_bad"},
            batch_key="batch",
            celltype_key="celltype",
            include_unintegrated=False,
            use_scib="sklearn",
        )
        assert isinstance(df, pd.DataFrame)
        assert "method" in df.columns
        assert "overall_score" in df.columns
        assert len(df) == 2
        # Should be sorted descending by overall_score
        assert df["overall_score"].iloc[0] >= df["overall_score"].iloc[1]

    def test_include_unintegrated(self):
        adata = _make_batched_adata()
        df = compare_integrations(
            adata,
            {"good": "X_good"},
            batch_key="batch",
            celltype_key="celltype",
            include_unintegrated=True,
            use_scib="sklearn",
        )
        assert any(df["method"] == "Unintegrated")

    def test_missing_embedding_skipped(self):
        adata = _make_batched_adata()
        df = compare_integrations(
            adata,
            {"good": "X_good", "missing": "X_nonexistent"},
            batch_key="batch",
            celltype_key="celltype",
            include_unintegrated=False,
            use_scib="sklearn",
        )
        assert len(df) == 1
        assert df["method"].iloc[0] == "good"

    def test_batch_bio_columns_present(self):
        adata = _make_batched_adata()
        df = compare_integrations(
            adata,
            {"good": "X_good"},
            batch_key="batch",
            celltype_key="celltype",
            include_unintegrated=False,
            use_scib="sklearn",
        )
        assert "batch_score" in df.columns
        assert "bio_score" in df.columns
        assert "overall_score" in df.columns

    def test_custom_weights(self):
        adata = _make_batched_adata()
        df1 = compare_integrations(
            adata,
            {"good": "X_good"},
            batch_key="batch",
            celltype_key="celltype",
            batch_weight=0.9,
            bio_weight=0.1,
            include_unintegrated=False,
            use_scib="sklearn",
        )
        df2 = compare_integrations(
            adata,
            {"good": "X_good"},
            batch_key="batch",
            celltype_key="celltype",
            batch_weight=0.1,
            bio_weight=0.9,
            include_unintegrated=False,
            use_scib="sklearn",
        )
        # Different weights should give different overall scores
        assert df1["overall_score"].iloc[0] != pytest.approx(
            df2["overall_score"].iloc[0], abs=0.001
        )


# ---------------------------------------------------------------------------
# sklearn fallback
# ---------------------------------------------------------------------------


class TestSklearnFallback:
    def test_explicit_sklearn(self):
        adata = _make_batched_adata()
        metrics = compute_integration_metrics(
            adata, "X_good", "batch", "celltype", use_scib="sklearn"
        )
        # Should have core sklearn metrics
        assert "asw_batch" in metrics
        assert "asw_celltype" in metrics
        assert "ari" in metrics
        assert "nmi" in metrics
        assert "pcr" in metrics

    def test_scib_required_without_install(self):
        # If scib-metrics is not installed, use_scib="scib" should raise
        try:
            import scib_metrics  # noqa: F401

            pytest.skip("scib-metrics is installed, cannot test missing case")
        except ImportError:
            adata = _make_batched_adata()
            with pytest.raises(ImportError, match="scib-metrics"):
                compute_integration_metrics(adata, "X_good", "batch", "celltype", use_scib="scib")


# ---------------------------------------------------------------------------
# Optional celltype_key
# ---------------------------------------------------------------------------


class TestOptionalCelltypeKey:
    def test_compute_integration_metrics_no_celltype(self):
        """Batch-only metrics returned when celltype_key is None."""
        adata = _make_batched_adata()
        metrics = compute_integration_metrics(
            adata, "X_good", "batch", celltype_key=None, use_scib="sklearn"
        )
        assert isinstance(metrics, dict)
        assert "asw_batch" in metrics
        assert "pcr" in metrics
        # Bio metrics should NOT be present
        assert "asw_celltype" not in metrics
        assert "ari" not in metrics
        assert "nmi" not in metrics

    def test_compare_integrations_no_celltype(self):
        """compare_integrations works with celltype_key=None."""
        adata = _make_batched_adata()
        df = compare_integrations(
            adata,
            {"good": "X_good"},
            batch_key="batch",
            celltype_key=None,
            include_unintegrated=False,
            use_scib="sklearn",
        )
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 1
        assert "overall_score" in df.columns
        assert "batch_score" in df.columns


# ---------------------------------------------------------------------------
# New integration methods
# ---------------------------------------------------------------------------


class TestRunCombat:
    def test_stores_pca_combat(self):
        """run_combat should store X_pca_combat in obsm."""
        import scanpy as sc

        # Use non-negative count data (randn has negatives → NaN after log1p)
        rng = np.random.RandomState(42)
        n_obs, n_vars, n_batches = 100, 50, 3
        X = rng.negative_binomial(5, 0.3, (n_obs, n_vars)).astype(np.float32)
        adata = AnnData(X=X)
        adata.obs["batch"] = [f"batch_{i % n_batches}" for i in range(n_obs)]

        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        from sc_tools.pp.integrate import run_combat

        run_combat(adata, batch_key="batch")
        assert "X_pca_combat" in adata.obsm
        assert adata.obsm["X_pca_combat"].shape[0] == adata.n_obs


class TestSoftDepHandling:
    def test_bbknn_import_error(self):
        """run_bbknn should raise ImportError with install hint if bbknn not available."""
        try:
            import bbknn  # noqa: F401

            pytest.skip("bbknn is installed")
        except ImportError:
            from sc_tools.pp.integrate import run_bbknn

            adata = _make_batched_adata()
            with pytest.raises(ImportError, match="bbknn"):
                run_bbknn(adata, batch_key="batch")

    def test_scanorama_import_error(self):
        """run_scanorama should raise ImportError with install hint if not available."""
        try:
            import scanorama  # noqa: F401

            pytest.skip("scanorama is installed")
        except ImportError:
            from sc_tools.pp.integrate import run_scanorama

            adata = _make_batched_adata()
            with pytest.raises(ImportError, match="scanorama"):
                run_scanorama(adata, batch_key="batch")


class TestRunIntegrationBenchmark:
    def test_returns_comparison_df(self):
        """run_integration_benchmark should return adata + comparison df."""
        import scanpy as sc

        from sc_tools.bm.integration import run_integration_benchmark

        # Use non-negative count data for combat compatibility
        rng = np.random.RandomState(42)
        n_obs, n_vars, n_batches, n_celltypes = 100, 50, 3, 4
        X = rng.negative_binomial(5, 0.3, (n_obs, n_vars)).astype(np.float32)
        adata = AnnData(X=X)
        adata.obs["batch"] = [f"batch_{i % n_batches}" for i in range(n_obs)]
        adata.obs["celltype"] = [f"type_{i % n_celltypes}" for i in range(n_obs)]
        adata.obsm["X_pca"] = rng.randn(n_obs, 10).astype(np.float32)

        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        # Only run methods that have no external deps
        adata_out, df = run_integration_benchmark(
            adata,
            modality="visium",
            batch_key="batch",
            methods=["combat", "pca"],
            use_scib="sklearn",
        )
        assert isinstance(df, pd.DataFrame)
        assert len(df) >= 2
        assert "overall_score" in df.columns
        assert "X_pca_combat" in adata_out.obsm


class TestExpandedEmbeddings:
    def test_known_embeddings_expanded(self):
        """_KNOWN_EMBEDDINGS should include new methods."""
        from sc_tools.qc.report_utils import _KNOWN_EMBEDDINGS

        assert "BBKNN" in _KNOWN_EMBEDDINGS
        assert "ComBat" in _KNOWN_EMBEDDINGS
        assert "Scanorama" in _KNOWN_EMBEDDINGS
        assert "scANVI" in _KNOWN_EMBEDDINGS
        assert "DestVI" in _KNOWN_EMBEDDINGS


# ---------------------------------------------------------------------------
# run_full_integration_workflow tests
# ---------------------------------------------------------------------------


class TestStratifiedSubsample:
    def test_subsample_preserves_batches(self):
        """Stratified subsample should include cells from all batches."""
        from sc_tools.bm.integration import _stratified_subsample

        adata = _make_batched_adata(n_obs=200)
        sub = _stratified_subsample(adata, "batch", n=50)
        assert sub.n_obs <= 50
        assert set(sub.obs["batch"].unique()) == set(adata.obs["batch"].unique())

    def test_subsample_noop_when_small(self):
        """Subsample returns copy when n >= n_obs."""
        from sc_tools.bm.integration import _stratified_subsample

        adata = _make_batched_adata(n_obs=40)
        sub = _stratified_subsample(adata, "batch", n=100)
        assert sub.n_obs == adata.n_obs


class TestResolveBestMethod:
    def test_uses_batch_score(self):
        """Should pick method with highest batch_score."""
        from sc_tools.bm.integration import _resolve_best_method

        df = pd.DataFrame(
            {
                "method": ["A", "B", "C"],
                "batch_score": [0.5, 0.8, 0.3],
                "overall_score": [0.9, 0.4, 0.7],
            }
        )
        assert _resolve_best_method(df) == "B"

    def test_falls_back_to_overall(self):
        """Falls back to overall_score when batch_score absent."""
        from sc_tools.bm.integration import _resolve_best_method

        df = pd.DataFrame(
            {
                "method": ["A", "B"],
                "overall_score": [0.3, 0.8],
            }
        )
        assert _resolve_best_method(df) == "B"


class TestRunFullIntegrationWorkflow:
    def test_basic_workflow(self, tmp_path):
        """Full workflow returns adata, comparison_df, and method name."""
        import scanpy as sc

        from sc_tools.bm.integration import run_full_integration_workflow

        rng = np.random.RandomState(42)
        n_obs, n_vars = 100, 50
        X = rng.negative_binomial(5, 0.3, (n_obs, n_vars)).astype(np.float32)
        adata = AnnData(X=X)
        adata.obs["batch"] = [f"batch_{i % 3}" for i in range(n_obs)]
        adata.obsm["X_pca"] = rng.randn(n_obs, 10).astype(np.float32)

        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        result_adata, df, best = run_full_integration_workflow(
            adata,
            modality="visium",
            batch_key="batch",
            methods=["combat", "pca"],
            output_dir=str(tmp_path),
            use_scib="sklearn",
            save_intermediates=True,
        )

        assert isinstance(df, pd.DataFrame)
        assert len(df) >= 2
        assert isinstance(best, str)

        # Check intermediates saved
        test_dir = tmp_path / "tmp" / "integration_test"
        assert test_dir.exists()
        assert len(list(test_dir.glob("*.h5ad"))) >= 1

        # Check method recorded
        method_file = tmp_path / "integration_method.txt"
        assert method_file.exists()
        assert method_file.read_text().strip() == best

    def test_no_save_intermediates(self, tmp_path):
        """save_intermediates=False skips file writes."""
        import scanpy as sc

        from sc_tools.bm.integration import run_full_integration_workflow

        rng = np.random.RandomState(42)
        n_obs, n_vars = 80, 40
        X = rng.negative_binomial(5, 0.3, (n_obs, n_vars)).astype(np.float32)
        adata = AnnData(X=X)
        adata.obs["batch"] = [f"batch_{i % 2}" for i in range(n_obs)]
        adata.obsm["X_pca"] = rng.randn(n_obs, 10).astype(np.float32)

        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        _, _, _ = run_full_integration_workflow(
            adata,
            batch_key="batch",
            methods=["pca"],
            output_dir=str(tmp_path),
            use_scib="sklearn",
            save_intermediates=False,
        )

        test_dir = tmp_path / "tmp" / "integration_test"
        assert not test_dir.exists()
