"""Unit tests for sc_tools.pp.projection -- subsample + kNN projection."""

from __future__ import annotations

import numpy as np
import pytest
from anndata import AnnData


@pytest.fixture()
def adata_clustered():
    """Synthetic AnnData with well-separated clusters and PCA."""
    np.random.seed(42)
    n_per_cluster = 200
    n_clusters = 3
    n_obs = n_per_cluster * n_clusters

    # Well-separated clusters in 10D PCA space
    pca = np.zeros((n_obs, 10), dtype="float32")
    labels = []
    for i in range(n_clusters):
        start = i * n_per_cluster
        end = (i + 1) * n_per_cluster
        pca[start:end] = np.random.randn(n_per_cluster, 10).astype("float32") + i * 10
        labels.extend([str(i)] * n_per_cluster)

    X = np.random.randn(n_obs, 50).astype("float32")
    adata = AnnData(X)
    adata.obsm["X_pca"] = pca
    adata.obs["leiden"] = labels
    adata.obs["leiden"] = adata.obs["leiden"].astype("category")
    adata.obs["library_id"] = (["L1"] * (n_obs // 2)) + (["L2"] * (n_obs // 2))
    adata.obs["library_id"] = adata.obs["library_id"].astype("category")
    # Fake UMAP
    adata.obsm["X_umap"] = pca[:, :2] + np.random.randn(n_obs, 2).astype("float32") * 0.1
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    return adata


# ---------------------------------------------------------------------------
# subsample_stratified
# ---------------------------------------------------------------------------


class TestSubsampleStratified:
    def test_proportions_within_tolerance(self, adata_clustered):
        from sc_tools.pp.projection import subsample_stratified

        idx = subsample_stratified(adata_clustered, n=300, stratify_key="library_id")
        assert len(idx) == 300
        # Check proportions roughly match (50/50 split)
        sub_obs = adata_clustered.obs.iloc[idx]
        counts = sub_obs["library_id"].value_counts()
        assert abs(counts["L1"] - 150) <= 10  # tolerance for rounding
        assert abs(counts["L2"] - 150) <= 10

    def test_deterministic_seeds(self, adata_clustered):
        from sc_tools.pp.projection import subsample_stratified

        idx1 = subsample_stratified(adata_clustered, n=200, random_state=42)
        idx2 = subsample_stratified(adata_clustered, n=200, random_state=42)
        np.testing.assert_array_equal(idx1, idx2)

    def test_n_geq_n_obs_returns_all(self, adata_clustered):
        from sc_tools.pp.projection import subsample_stratified

        idx = subsample_stratified(adata_clustered, n=10000)
        assert len(idx) == adata_clustered.n_obs

    def test_no_stratify_key_random(self, adata_clustered):
        from sc_tools.pp.projection import subsample_stratified

        idx = subsample_stratified(adata_clustered, n=100, stratify_key=None)
        assert len(idx) == 100
        assert len(np.unique(idx)) == 100  # no duplicates

    def test_zero_count_groups_skipped(self):
        """Groups with 0 cells in a categorical should be skipped."""
        import pandas as pd

        from sc_tools.pp.projection import subsample_stratified

        adata = AnnData(np.random.randn(100, 5).astype("float32"))
        # Category "C" exists but has 0 cells
        adata.obs["batch"] = pd.Categorical(["A"] * 60 + ["B"] * 40, categories=["A", "B", "C"])
        idx = subsample_stratified(adata, n=50, stratify_key="batch")
        assert len(idx) == 50

    def test_single_sample_fallback(self):
        """Single-sample dataset should work without error."""
        from sc_tools.pp.projection import subsample_stratified

        adata = AnnData(np.random.randn(200, 5).astype("float32"))
        adata.obs["batch"] = "single"
        idx = subsample_stratified(adata, n=50, stratify_key="batch")
        assert len(idx) == 50


# ---------------------------------------------------------------------------
# build_knn_index
# ---------------------------------------------------------------------------


class TestBuildKnnIndex:
    def test_builds_index(self, adata_clustered):
        from sc_tools.pp.projection import build_knn_index

        # Should build a BallTree index (FAISS unlikely in test env)
        index = build_knn_index(adata_clustered, use_rep="X_pca")
        assert index is not None

    def test_faiss_missing_message(self, adata_clustered, monkeypatch):
        """Should log info when FAISS unavailable, not error."""
        import builtins

        from sc_tools.pp.projection import build_knn_index

        _real_import = builtins.__import__

        def mock_import(name, *args, **kwargs):
            if "faiss" in name:
                raise ImportError("mocked")
            return _real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", mock_import)
        index = build_knn_index(adata_clustered, use_rep="X_pca")
        assert index is not None  # fell back to BallTree


# ---------------------------------------------------------------------------
# project_labels
# ---------------------------------------------------------------------------


class TestProjectLabels:
    def test_well_separated_clusters(self, adata_clustered):
        """Well-separated clusters should have ~100% confidence."""
        from sc_tools.pp.projection import (
            build_knn_index,
            project_labels,
            subsample_stratified,
        )

        idx = subsample_stratified(adata_clustered, n=300, random_state=42)
        adata_sub = adata_clustered[idx].copy()
        index = build_knn_index(adata_sub, use_rep="X_pca")

        project_labels(adata_clustered, index, adata_sub, label_key="leiden", k=15)

        assert "leiden_projected" in adata_clustered.obs.columns
        assert "leiden_confidence" in adata_clustered.obs.columns
        # Confidence should be in [0, 1]
        conf = adata_clustered.obs["leiden_confidence"].values
        assert np.all(conf >= 0) and np.all(conf <= 1)
        # Well-separated -> high confidence
        assert np.median(conf) > 0.8


# ---------------------------------------------------------------------------
# project_umap
# ---------------------------------------------------------------------------


class TestProjectUmap:
    def test_umap_projection(self, adata_clustered):
        from sc_tools.pp.projection import (
            build_knn_index,
            project_umap,
            subsample_stratified,
        )

        idx = subsample_stratified(adata_clustered, n=300, random_state=42)
        adata_sub = adata_clustered[idx].copy()
        index = build_knn_index(adata_sub, use_rep="X_pca")

        project_umap(adata_clustered, index, adata_sub, k=15)

        assert "X_umap" in adata_clustered.obsm
        assert adata_clustered.obsm["X_umap"].shape == (600, 2)


# ---------------------------------------------------------------------------
# project_representation
# ---------------------------------------------------------------------------


class TestProjectRepresentation:
    def test_correct_shape(self, adata_clustered):
        from sc_tools.pp.projection import (
            build_knn_index,
            project_representation,
            subsample_stratified,
        )

        idx = subsample_stratified(adata_clustered, n=300, random_state=42)
        adata_sub = adata_clustered[idx].copy()
        # Fake a corrected representation
        adata_sub.obsm["X_pca_harmony"] = adata_sub.obsm["X_pca"] + 0.1
        index = build_knn_index(adata_sub, use_rep="X_pca")

        project_representation(adata_clustered, index, adata_sub, rep_key="X_pca_harmony", k=15)

        assert "X_pca_harmony" in adata_clustered.obsm
        assert adata_clustered.obsm["X_pca_harmony"].shape == (600, 10)


# ---------------------------------------------------------------------------
# _store_projection_info
# ---------------------------------------------------------------------------


class TestStoreProjectionInfo:
    def test_required_keys(self, adata_clustered):
        from sc_tools.pp.projection import _store_projection_info
        from sc_tools.pp.strategy import SubsampleContext

        ctx = SubsampleContext(
            subsample_idx=np.array([0, 1, 2]),
            subsample_n=3,
            projection_k=15,
            random_state=42,
            stratify_key="library_id",
            use_rep="X_pca",
        )

        _store_projection_info(adata_clustered, ctx, strategy_name="large")

        assert "_projection_info" in adata_clustered.uns
        info = adata_clustered.uns["_projection_info"]
        assert info["strategy"] == "large"
        assert info["subsample_n"] == 3
        assert info["projection_k"] == 15
        assert info["random_state"] == 42
        assert info["use_rep"] == "X_pca"

        # obs mask
        assert "_sc_tools_is_representative" in adata_clustered.obs.columns
        mask = adata_clustered.obs["_sc_tools_is_representative"]
        assert mask.sum() == 3
