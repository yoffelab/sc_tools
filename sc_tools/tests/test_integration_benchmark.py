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

    def test_missing_celltype_key_raises(self):
        adata = _make_batched_adata()
        with pytest.raises(KeyError, match="missing_ct"):
            compute_integration_metrics(adata, "X_good", "batch", "missing_ct", use_scib="sklearn")

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
