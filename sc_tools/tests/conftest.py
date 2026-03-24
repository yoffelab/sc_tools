"""Shared test fixtures for sc_tools tests."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData


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
