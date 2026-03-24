"""Tests for sc_tools.io.estimate -- memory estimation from h5 metadata."""

from __future__ import annotations

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
    import pandas as pd

    n_obs, n_vars = 50, 20
    X = np.random.default_rng(42).standard_normal((n_obs, n_vars)).astype(np.float32)
    obs_df = pd.DataFrame(
        {"batch": [f"b{i % 3}" for i in range(n_obs)]},
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    var_df = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])
    adata = ad.AnnData(X=X, obs=obs_df, var=var_df)
    path = tmp_path / "dense.h5ad"
    adata.write_h5ad(path)
    return path


@pytest.fixture()
def tmp_sparse_h5ad(tmp_path):
    """Create a small sparse CSR h5ad for testing."""
    import pandas as pd

    n_obs, n_vars = 50, 20
    rng = np.random.default_rng(42)
    dense = rng.standard_normal((n_obs, n_vars)).astype(np.float32)
    dense[dense < 0.5] = 0
    X = sp.csr_matrix(dense)
    obs_df = pd.DataFrame(
        {"batch": [f"b{i % 3}" for i in range(n_obs)]},
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    var_df = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])
    adata = ad.AnnData(X=X, obs=obs_df, var=var_df)
    path = tmp_path / "sparse.h5ad"
    adata.write_h5ad(path)
    return path


@pytest.fixture()
def tmp_zero_genes_h5ad(tmp_path):
    """Create an h5ad with 0 genes."""
    import pandas as pd

    n_obs, n_vars = 10, 0
    X = np.empty((n_obs, n_vars), dtype=np.float32)
    obs_df = pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)])
    var_df = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])
    adata = ad.AnnData(X=X, obs=obs_df, var=var_df)
    path = tmp_path / "zero_genes.h5ad"
    adata.write_h5ad(path)
    return path


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestEstimateFromH5:
    def test_estimate_dense(self, tmp_dense_h5ad):
        from sc_tools.io.estimate import estimate_from_h5

        result = estimate_from_h5(tmp_dense_h5ad)
        assert isinstance(result, dict)
        assert "x_mb" in result
        assert "obs_mb" in result
        assert "base_mb" in result
        assert "estimated_peak_mb" in result
        # Peak should be higher than base due to overhead multiplier
        assert result["estimated_peak_mb"] > result["base_mb"]

    def test_estimate_sparse(self, tmp_dense_h5ad, tmp_sparse_h5ad):
        from sc_tools.io.estimate import estimate_from_h5

        dense_est = estimate_from_h5(tmp_dense_h5ad)
        sparse_est = estimate_from_h5(tmp_sparse_h5ad)
        # Sparse x_mb should be less than dense equivalent
        assert sparse_est["x_mb"] < dense_est["x_mb"]

    def test_estimate_keys(self, tmp_dense_h5ad):
        from sc_tools.io.estimate import estimate_from_h5

        result = estimate_from_h5(tmp_dense_h5ad)
        expected_keys = {
            "n_obs",
            "n_vars",
            "x_dtype",
            "x_sparse",
            "x_mb",
            "obs_mb",
            "var_mb",
            "base_mb",
            "estimated_peak_mb",
        }
        assert expected_keys.issubset(set(result.keys()))

    def test_estimate_zero_genes(self, tmp_zero_genes_h5ad):
        from sc_tools.io.estimate import estimate_from_h5

        result = estimate_from_h5(tmp_zero_genes_h5ad)
        # With zero genes, estimated_peak_mb should be 0 or small positive
        assert result["estimated_peak_mb"] >= 0
        assert result["x_mb"] == 0.0
