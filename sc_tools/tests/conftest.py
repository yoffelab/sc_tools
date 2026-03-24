"""Shared test fixtures for sc_tools tests."""

from __future__ import annotations

import numpy as np
import pandas as pd
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


@pytest.fixture
def adata_multi_subject() -> AnnData:
    """6 subjects (3 control, 3 treatment), 50 cells/subject, 200 genes, 3 celltypes.

    Batch 0 = subjects 0-2, batch 1 = subjects 3-5.
    NOT confounded: each batch has both conditions (subject 0,1 = control batch 0;
    subject 2 = treatment batch 0; subject 3 = control batch 1; subject 4,5 = treatment batch 1).
    """
    rng = np.random.default_rng(42)
    n_subjects = 6
    n_cells_per_subject = 50
    n_genes = 200
    n_obs = n_subjects * n_cells_per_subject

    X = rng.negative_binomial(5, 0.3, (n_obs, n_genes)).astype("float32")

    subjects = np.repeat(
        [f"subject_{i}" for i in range(n_subjects)], n_cells_per_subject
    )
    # Conditions: subjects 0,1,3 = control; subjects 2,4,5 = treatment
    condition_map = {
        "subject_0": "control",
        "subject_1": "control",
        "subject_2": "treatment",
        "subject_3": "control",
        "subject_4": "treatment",
        "subject_5": "treatment",
    }
    conditions = [condition_map[s] for s in subjects]

    # Batches: subjects 0-2 = batch_0, subjects 3-5 = batch_1
    batch_map = {
        "subject_0": "batch_0",
        "subject_1": "batch_0",
        "subject_2": "batch_0",
        "subject_3": "batch_1",
        "subject_4": "batch_1",
        "subject_5": "batch_1",
    }
    batches = [batch_map[s] for s in subjects]

    libraries = [f"lib_{s}" for s in subjects]
    celltypes = rng.choice(
        [f"type_{i}" for i in range(3)], n_obs
    )

    obs = pd.DataFrame(
        {
            "subject_id": pd.Categorical(subjects),
            "library_id": pd.Categorical(libraries),
            "condition": pd.Categorical(conditions),
            "batch": pd.Categorical(batches),
            "celltype": pd.Categorical(celltypes),
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )

    return AnnData(
        X,
        obs=obs,
        var=pd.DataFrame(index=[f"Gene{i}" for i in range(n_genes)]),
    )


@pytest.fixture
def adata_panel() -> AnnData:
    """Panel-sized AnnData: 6 subjects, 50 cells/subject, 40 genes (IMC-like).

    Same subject/condition/batch structure as adata_multi_subject but n_vars=40.
    """
    rng = np.random.default_rng(42)
    n_subjects = 6
    n_cells_per_subject = 50
    n_genes = 40
    n_obs = n_subjects * n_cells_per_subject

    X = rng.negative_binomial(5, 0.3, (n_obs, n_genes)).astype("float32")

    subjects = np.repeat(
        [f"subject_{i}" for i in range(n_subjects)], n_cells_per_subject
    )
    libraries = [f"lib_{s}" for s in subjects]

    obs = pd.DataFrame(
        {
            "subject_id": pd.Categorical(subjects),
            "library_id": pd.Categorical(libraries),
            "condition": pd.Categorical(
                ["control"] * (n_obs // 2) + ["treatment"] * (n_obs // 2)
            ),
            "batch": pd.Categorical(
                ["batch_0"] * (n_obs // 2) + ["batch_1"] * (n_obs // 2)
            ),
            "celltype": pd.Categorical(
                rng.choice([f"type_{i}" for i in range(3)], n_obs)
            ),
            "leiden": pd.Categorical([str(i % 4) for i in range(n_obs)]),
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )

    return AnnData(
        X,
        obs=obs,
        var=pd.DataFrame(index=[f"Marker{i}" for i in range(n_genes)]),
    )
