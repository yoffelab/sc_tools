"""
Unit tests for sc_tools.tl.deconvolution.

Tests use synthetic data and mock backends -- no real tangram/cell2location
needed.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from sc_tools.tl.deconvolution import (
    _BACKENDS,
    deconvolution,
    extract_reference_profiles,
    get_backend,
    register_backend,
    select_signature_genes,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_sc_adata(n_obs=200, n_vars=100, n_celltypes=4, seed=42):
    """Synthetic scRNA-seq reference with cell-type labels and batch."""
    rng = np.random.default_rng(seed)
    X = rng.negative_binomial(5, 0.3, (n_obs, n_vars)).astype(np.float32)
    var_names = [f"Gene{i}" for i in range(n_vars)]
    obs = pd.DataFrame(
        {
            "celltype": rng.choice([f"CT{i}" for i in range(n_celltypes)], n_obs),
            "batch": rng.choice(["B1", "B2"], n_obs),
        },
        index=[f"sc_cell_{i}" for i in range(n_obs)],
    )
    adata = sc.AnnData(X, obs=obs, var=pd.DataFrame(index=var_names))
    adata.raw = adata.copy()
    return adata


def _make_spatial_adata(n_obs=300, n_vars=100, n_libs=3, seed=99):
    """Synthetic spatial AnnData with library_id."""
    rng = np.random.default_rng(seed)
    X = rng.negative_binomial(5, 0.3, (n_obs, n_vars)).astype(np.float32)
    var_names = [f"Gene{i}" for i in range(n_vars)]
    obs = pd.DataFrame(
        {"library_id": rng.choice([f"lib{i}" for i in range(n_libs)], n_obs)},
        index=[f"spot_{i}" for i in range(n_obs)],
    )
    adata = sc.AnnData(X, obs=obs, var=pd.DataFrame(index=var_names))
    adata.raw = adata.copy()
    return adata


# ---------------------------------------------------------------------------
# Mock backend for testing orchestration
# ---------------------------------------------------------------------------

class MockBackend:
    """Backend that returns uniform proportions without any real model."""

    @staticmethod
    def run(
        sc_adata,
        spatial_adata_lib,
        shared_genes,
        celltype_key,
        *,
        use_gpu=False,
        reference_profiles=None,
        logger_instance=None,
        **kwargs,
    ):
        n_spots = spatial_adata_lib.n_obs
        n_cts = kwargs.get("n_celltypes", 4)
        rng = np.random.default_rng(0)
        props = rng.dirichlet(np.ones(n_cts), size=n_spots).astype(np.float32)
        return props


# ---------------------------------------------------------------------------
# Tests: extract_reference_profiles
# ---------------------------------------------------------------------------

class TestExtractReferenceProfiles:
    def test_basic(self):
        sc_adata = _make_sc_adata(n_obs=100, n_vars=50, n_celltypes=3)
        df = extract_reference_profiles(sc_adata, "celltype")
        assert isinstance(df, pd.DataFrame)
        assert df.shape == (50, 3)
        assert set(df.columns) == {"CT0", "CT1", "CT2"}
        assert list(df.index) == [f"Gene{i}" for i in range(50)]

    def test_gene_subset(self):
        sc_adata = _make_sc_adata(n_obs=100, n_vars=50, n_celltypes=3)
        genes = ["Gene0", "Gene1", "Gene2", "NonExistent"]
        df = extract_reference_profiles(sc_adata, "celltype", genes=genes)
        assert df.shape[0] == 3  # only existing genes
        assert "NonExistent" not in df.index

    def test_qc_labels(self):
        sc_adata = _make_sc_adata(n_obs=100, n_vars=50, n_celltypes=4)
        df_full = extract_reference_profiles(sc_adata, "celltype")
        df_filtered = extract_reference_profiles(
            sc_adata, "celltype", qc_labels=["CT0"]
        )
        assert "CT0" not in df_filtered.columns
        assert df_filtered.shape[1] == df_full.shape[1] - 1

    def test_cache(self, tmp_path):
        sc_adata = _make_sc_adata(n_obs=100, n_vars=50, n_celltypes=3)
        cache = tmp_path / "ref_profiles.pkl"
        df1 = extract_reference_profiles(sc_adata, "celltype", cache_path=cache)
        assert cache.exists()
        df2 = extract_reference_profiles(sc_adata, "celltype", cache_path=cache)
        pd.testing.assert_frame_equal(df1, df2)


# ---------------------------------------------------------------------------
# Tests: backend registry
# ---------------------------------------------------------------------------

class TestBackendRegistry:
    def test_builtin_backends(self):
        assert "tangram" in _BACKENDS
        assert "cell2location" in _BACKENDS
        assert "destvi" in _BACKENDS

    def test_register_and_get(self):
        register_backend("mock_test", MockBackend)
        assert get_backend("mock_test") is MockBackend
        # Cleanup
        del _BACKENDS["mock_test"]

    def test_unknown_backend_raises(self):
        with pytest.raises(ValueError, match="Unknown deconvolution method"):
            get_backend("nonexistent_method")


# ---------------------------------------------------------------------------
# Tests: deconvolution() with mock backend
# ---------------------------------------------------------------------------

class TestDeconvolutionOrchestrator:
    def setup_method(self):
        register_backend("mock", MockBackend)

    def teardown_method(self):
        if "mock" in _BACKENDS:
            del _BACKENDS["mock"]

    def test_basic_deconvolution(self):
        sc_adata = _make_sc_adata(n_obs=80, n_vars=100, n_celltypes=4)
        spatial = _make_spatial_adata(n_obs=150, n_vars=100, n_libs=2)

        result = deconvolution(
            spatial,
            sc_adata,
            method="mock",
            celltype_key="celltype",
            spatial_batch_key="library_id",
            sc_batch_key="batch",
            method_kwargs={"n_celltypes": 4},
        )

        assert "cell_type_proportions" in result.obsm
        props = result.obsm["cell_type_proportions"]
        assert isinstance(props, pd.DataFrame)
        assert props.shape[0] == 150
        assert props.shape[1] == 4
        assert "mock_argmax" in result.obs.columns

    def test_argmax_labels(self):
        sc_adata = _make_sc_adata(n_obs=50, n_vars=100, n_celltypes=3)
        spatial = _make_spatial_adata(n_obs=60, n_vars=100, n_libs=1)

        result = deconvolution(
            spatial,
            sc_adata,
            method="mock",
            celltype_key="celltype",
            spatial_batch_key="library_id",
            method_kwargs={"n_celltypes": 3},
        )

        argmax_col = result.obs["mock_argmax"]
        assert len(argmax_col) == 60
        # All values should be valid cell type names or generic labels
        assert all(isinstance(v, str) for v in argmax_col)

    def test_no_batch_key(self):
        """When spatial_batch_key is missing, fall back to single-batch mode."""
        sc_adata = _make_sc_adata(n_obs=50, n_vars=100, n_celltypes=3)
        spatial = _make_spatial_adata(n_obs=30, n_vars=100, n_libs=1)
        del spatial.obs["library_id"]

        result = deconvolution(
            spatial,
            sc_adata,
            method="mock",
            celltype_key="celltype",
            spatial_batch_key="library_id",
            method_kwargs={"n_celltypes": 3},
        )
        assert "cell_type_proportions" in result.obsm

    def test_output_file(self, tmp_path):
        sc_adata = _make_sc_adata(n_obs=50, n_vars=100, n_celltypes=3)
        spatial = _make_spatial_adata(n_obs=30, n_vars=100, n_libs=1)
        out = tmp_path / "result.h5ad"

        deconvolution(
            spatial,
            sc_adata,
            method="mock",
            celltype_key="celltype",
            spatial_batch_key="library_id",
            method_kwargs={"n_celltypes": 3},
            output_file=out,
        )
        assert out.exists()

    def test_proportion_alignment(self):
        """Proportions must be aligned to the original obs order."""
        sc_adata = _make_sc_adata(n_obs=50, n_vars=100, n_celltypes=2)
        spatial = _make_spatial_adata(n_obs=100, n_vars=100, n_libs=3)

        result = deconvolution(
            spatial,
            sc_adata,
            method="mock",
            celltype_key="celltype",
            spatial_batch_key="library_id",
            method_kwargs={"n_celltypes": 2},
        )

        props = result.obsm["cell_type_proportions"]
        assert list(props.index) == list(spatial.obs_names)


# ---------------------------------------------------------------------------
# Tests: select_signature_genes
# ---------------------------------------------------------------------------

class TestSelectSignatureGenes:
    def test_basic(self):
        sc_adata = _make_sc_adata(n_obs=100, n_vars=100, n_celltypes=3)
        genes = select_signature_genes(
            sc_adata, "celltype", "batch", n_genes_max=50, skip_hvg=True
        )
        assert isinstance(genes, list)
        assert len(genes) > 0
        assert len(genes) <= 100  # cannot exceed n_vars
        assert all(g in sc_adata.var_names for g in genes)

    def test_cache(self, tmp_path):
        sc_adata = _make_sc_adata(n_obs=100, n_vars=100, n_celltypes=3)
        genes1 = select_signature_genes(
            sc_adata,
            "celltype",
            "batch",
            n_genes_max=50,
            skip_hvg=True,
            cache_dir=str(tmp_path),
            sc_data_file="fake_path.h5ad",
        )
        # Second call should load from cache
        genes2 = select_signature_genes(
            sc_adata,
            "celltype",
            "batch",
            n_genes_max=50,
            skip_hvg=True,
            cache_dir=str(tmp_path),
            sc_data_file="fake_path.h5ad",
        )
        assert genes1 == genes2
