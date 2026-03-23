"""Shared test fixtures for sc_tools CLI integration tests (TST-05)."""

from __future__ import annotations

import numpy as np
import pytest
from anndata import AnnData


@pytest.fixture
def adata_100() -> AnnData:
    """100-cell AnnData with batch, celltype, spatial, and embeddings for CLI integration tests.

    Has two libraries (L0, L1) with 50 cells each, two cell types, spatial coordinates,
    and a PCA-like embedding in obsm. Gene names include MT-prefixed genes for QC metrics.
    """
    rng = np.random.default_rng(42)
    n_obs, n_vars = 100, 200

    # Raw count matrix (negative binomial for realism)
    X = rng.negative_binomial(5, 0.3, (n_obs, n_vars)).astype("float32")

    # Gene names: 3 MT genes + 197 regular genes
    var_names = [f"MT-{i}" for i in range(3)] + [f"GENE{i}" for i in range(3, n_vars)]

    adata = AnnData(X)
    adata.var_names = var_names
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]

    # obs columns matching pipeline requirements
    adata.obs["library_id"] = [f"L{i % 2}" for i in range(n_obs)]
    adata.obs["sample"] = adata.obs["library_id"]
    adata.obs["raw_data_dir"] = [f"/data/L{i % 2}" for i in range(n_obs)]
    adata.obs["celltype"] = [f"type_{i % 2}" for i in range(n_obs)]
    adata.obs["celltype_broad"] = [f"broad_{i % 2}" for i in range(n_obs)]
    adata.obs["leiden"] = [str(i % 3) for i in range(n_obs)]

    # obsm: spatial coordinates + PCA-like embedding
    adata.obsm["spatial"] = rng.random((n_obs, 2)).astype("float32")
    adata.obsm["X_pca"] = rng.standard_normal((n_obs, 50)).astype("float32")
    adata.obsm["X_scvi"] = rng.standard_normal((n_obs, 10)).astype("float32")

    # uns: modality hint for auto-detection
    adata.uns["modality"] = "visium"

    return adata


@pytest.fixture
def adata_100_h5ad(adata_100: AnnData, tmp_path) -> str:
    """Write adata_100 to a temporary h5ad file, return the path as string."""
    path = tmp_path / "test_data.h5ad"
    adata_100.write_h5ad(path)
    return str(path)


@pytest.fixture
def adata_100_preprocess_checkpoint(adata_100: AnnData, tmp_path) -> str:
    """AnnData that passes 'preprocess' phase validation (has integration rep + clusters + raw backup)."""
    adata_100.raw = adata_100.copy()
    path = tmp_path / "test_preprocessed.h5ad"
    adata_100.write_h5ad(path)
    return str(path)
