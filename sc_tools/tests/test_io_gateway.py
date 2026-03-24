"""Tests for sc_tools.io gateway, metadata reader, and tiered loading."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def tmp_dense_h5ad(tmp_path):
    """Create a small dense h5ad for testing."""
    n_obs, n_vars = 50, 20
    X = np.random.default_rng(42).standard_normal((n_obs, n_vars)).astype(np.float32)
    obs = {
        "batch": [f"b{i % 3}" for i in range(n_obs)],
        "celltype": [f"ct{i % 5}" for i in range(n_obs)],
    }
    import pandas as pd

    obs_df = pd.DataFrame(obs, index=[f"cell_{i}" for i in range(n_obs)])
    var_df = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])
    adata = ad.AnnData(X=X, obs=obs_df, var=var_df)
    adata.obsm["X_pca"] = np.random.default_rng(0).standard_normal((n_obs, 10)).astype(
        np.float32
    )
    path = tmp_path / "dense.h5ad"
    adata.write_h5ad(path)
    return path


@pytest.fixture()
def tmp_sparse_h5ad(tmp_path):
    """Create a small sparse CSR h5ad for testing."""
    n_obs, n_vars = 50, 20
    rng = np.random.default_rng(42)
    dense = rng.standard_normal((n_obs, n_vars)).astype(np.float32)
    dense[dense < 0.5] = 0  # make it sparse-ish
    X = sp.csr_matrix(dense)
    import pandas as pd

    obs_df = pd.DataFrame(
        {
            "batch": [f"b{i % 3}" for i in range(n_obs)],
            "celltype": [f"ct{i % 5}" for i in range(n_obs)],
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    var_df = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])
    adata = ad.AnnData(X=X, obs=obs_df, var=var_df)
    adata.obsm["X_pca"] = rng.standard_normal((n_obs, 10)).astype(np.float32)
    path = tmp_path / "sparse.h5ad"
    adata.write_h5ad(path)
    return path


# ---------------------------------------------------------------------------
# Metadata reader tests
# ---------------------------------------------------------------------------


class TestReadH5adMetadata:
    def test_t1_metadata_read(self, tmp_dense_h5ad):
        from sc_tools.io.metadata import read_h5ad_metadata

        meta = read_h5ad_metadata(tmp_dense_h5ad)
        assert isinstance(meta, dict)
        assert meta["n_obs"] == 50
        assert meta["n_vars"] == 20
        assert meta["x_dtype"] == "float32"
        assert meta["x_sparse"] is False
        assert "batch" in meta["obs_columns"]
        assert "celltype" in meta["obs_columns"]
        assert "X_pca" in meta["obsm_keys"]

    def test_t1_metadata_sparse(self, tmp_sparse_h5ad):
        from sc_tools.io.metadata import read_h5ad_metadata

        meta = read_h5ad_metadata(tmp_sparse_h5ad)
        assert meta["x_sparse"] is True
        assert meta["x_dtype"] == "float32"
        assert meta["n_obs"] == 50
        assert meta["n_vars"] == 20
        assert "nnz" in meta


# ---------------------------------------------------------------------------
# IOGateway tests
# ---------------------------------------------------------------------------


class TestIOGateway:
    def test_t1_returns_dict_not_anndata(self, tmp_dense_h5ad):
        from sc_tools.io.gateway import DataTier, IOGateway

        result = IOGateway().read(tmp_dense_h5ad, DataTier.T1_METADATA)
        assert isinstance(result, dict)
        assert not isinstance(result, ad.AnnData)
        assert result["n_obs"] == 50

    def test_t2_backed_read(self, tmp_dense_h5ad):
        from sc_tools.io.gateway import DataTier, IOGateway

        result = IOGateway().read(tmp_dense_h5ad, DataTier.T2_SUMMARY)
        assert isinstance(result, ad.AnnData)
        assert result.isbacked

    def test_t3_full_read(self, tmp_dense_h5ad):
        from sc_tools.io.gateway import DataTier, IOGateway

        result = IOGateway().read(tmp_dense_h5ad, DataTier.T3_FULL)
        assert isinstance(result, ad.AnnData)
        assert not result.isbacked

    def test_memory_guard_blocks(self, tmp_dense_h5ad):
        from sc_tools.errors import SCToolsRuntimeError
        from sc_tools.io.gateway import DataTier, IOGateway

        # Mock estimate to return a large value, and psutil to return small available
        large_estimate = {
            "n_obs": 2_500_000,
            "n_vars": 30_000,
            "x_dtype": "float32",
            "x_sparse": False,
            "x_mb": 50_000.0,
            "obs_mb": 500.0,
            "var_mb": 3.0,
            "base_mb": 50_503.0,
            "estimated_peak_mb": 101_006.0,
        }
        mock_vmem = MagicMock()
        mock_vmem.available = 8 * 1024 * 1024 * 1024  # 8 GB

        with patch("sc_tools.io.gateway.estimate_from_h5", return_value=large_estimate):
            with patch("psutil.virtual_memory", return_value=mock_vmem):
                with pytest.raises(SCToolsRuntimeError, match="exceeds 80%"):
                    IOGateway().read(tmp_dense_h5ad, DataTier.T3_FULL)

    def test_force_override(self, tmp_dense_h5ad):
        from sc_tools.io.gateway import DataTier, IOGateway

        # Mock estimate to return a large value, and psutil to return small available
        large_estimate = {
            "estimated_peak_mb": 101_006.0,
        }
        mock_vmem = MagicMock()
        mock_vmem.available = 8 * 1024 * 1024 * 1024  # 8 GB

        with patch("sc_tools.io.gateway.estimate_from_h5", return_value=large_estimate):
            with patch("psutil.virtual_memory", return_value=mock_vmem):
                result = IOGateway().read(tmp_dense_h5ad, DataTier.T3_FULL, force=True)
                assert isinstance(result, ad.AnnData)
                assert not result.isbacked
