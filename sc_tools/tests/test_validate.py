"""Unit tests for sc_tools.validate checkpoint validation."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from sc_tools.validate import (
    CheckpointValidationError,
    validate_checkpoint,
    validate_file,
    validate_p1,
    validate_p2,
    validate_p3,
    validate_p4,
    validate_p35,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def adata_p1():
    """Minimal valid Phase 1 AnnData."""
    n_obs, n_vars = 50, 100
    X = np.random.poisson(5, (n_obs, n_vars)).astype(np.float32)
    obs = pd.DataFrame(
        {
            "sample": [f"sample_{i % 3}" for i in range(n_obs)],
            "raw_data_dir": [f"/data/sample_{i % 3}" for i in range(n_obs)],
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obsm["spatial"] = np.random.rand(n_obs, 2)
    return adata


@pytest.fixture
def adata_p2(adata_p1):
    """Minimal valid Phase 2 AnnData (p1 + clinical metadata)."""
    adata_p1.obs["diagnosis"] = "cancer"
    adata_p1.obs["age"] = 55
    return adata_p1


@pytest.fixture
def adata_p3():
    """Minimal valid Phase 3 AnnData."""
    n_obs, n_vars = 50, 100
    X = np.random.rand(n_obs, n_vars).astype(np.float32)
    obs = pd.DataFrame(
        {"leiden": [str(i % 5) for i in range(n_obs)]},
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obsm["X_scvi"] = np.random.rand(n_obs, 10)
    # Backup raw
    adata.raw = adata.copy()
    return adata


@pytest.fixture
def adata_p35(adata_p3):
    """Minimal valid Phase 3.5b AnnData."""
    n_obs = adata_p3.n_obs
    adata_p3.obsm["signature_score"] = pd.DataFrame(
        {"Hallmark/HYPOXIA": np.random.rand(n_obs)},
        index=adata_p3.obs_names,
    )
    adata_p3.obsm["signature_score_z"] = pd.DataFrame(
        {"Hallmark/HYPOXIA": np.random.randn(n_obs)},
        index=adata_p3.obs_names,
    )
    adata_p3.uns["signature_score_report"] = {"Hallmark/HYPOXIA": {"n_present": 20, "n_missing": 5}}
    return adata_p3


@pytest.fixture
def adata_p4(adata_p35):
    """Minimal valid Phase 4 AnnData."""
    adata_p35.obs["celltype"] = "Epithelial"
    adata_p35.obs["celltype_broad"] = "Epithelial"
    return adata_p35


# ---------------------------------------------------------------------------
# validate_p1
# ---------------------------------------------------------------------------


class TestValidateP1:
    def test_valid(self, adata_p1):
        issues = validate_p1(adata_p1)
        assert issues == []

    def test_missing_sample(self, adata_p1):
        del adata_p1.obs["sample"]
        issues = validate_p1(adata_p1)
        assert any("obs['sample']" in i for i in issues)

    def test_missing_raw_data_dir(self, adata_p1):
        del adata_p1.obs["raw_data_dir"]
        issues = validate_p1(adata_p1)
        assert any("obs['raw_data_dir']" in i for i in issues)

    def test_missing_spatial(self, adata_p1):
        del adata_p1.obsm["spatial"]
        issues = validate_p1(adata_p1)
        assert any("obsm['spatial']" in i for i in issues)

    def test_negative_X(self, adata_p1):
        adata_p1.X[0, 0] = -1.0
        issues = validate_p1(adata_p1)
        assert any("negative" in i for i in issues)

    def test_fix_batch_to_raw_data_dir(self, adata_p1):
        del adata_p1.obs["raw_data_dir"]
        adata_p1.obs["batch"] = "batch_val"
        issues = validate_p1(adata_p1, fix=True)
        assert issues == []
        assert "raw_data_dir" in adata_p1.obs.columns

    def test_fix_hint_without_fix_flag(self, adata_p1):
        del adata_p1.obs["raw_data_dir"]
        adata_p1.obs["batch"] = "batch_val"
        issues = validate_p1(adata_p1, fix=False)
        assert any("fix=True" in i for i in issues)


# ---------------------------------------------------------------------------
# validate_p2
# ---------------------------------------------------------------------------


class TestValidateP2:
    def test_valid(self, adata_p2):
        issues = validate_p2(adata_p2)
        assert issues == []

    def test_no_clinical_metadata(self, adata_p1):
        """p1-valid but no extra clinical columns -> p2 fails."""
        issues = validate_p2(adata_p1)
        assert any("clinical metadata" in i for i in issues)

    def test_inherits_p1_issues(self, adata_p2):
        del adata_p2.obs["sample"]
        issues = validate_p2(adata_p2)
        assert any("obs['sample']" in i for i in issues)


# ---------------------------------------------------------------------------
# validate_p3
# ---------------------------------------------------------------------------


class TestValidateP3:
    def test_valid(self, adata_p3):
        issues = validate_p3(adata_p3)
        assert issues == []

    def test_missing_integration_rep(self, adata_p3):
        del adata_p3.obsm["X_scvi"]
        issues = validate_p3(adata_p3)
        assert any("integration representation" in i for i in issues)

    def test_missing_cluster_col(self, adata_p3):
        del adata_p3.obs["leiden"]
        issues = validate_p3(adata_p3)
        assert any("cluster column" in i for i in issues)

    def test_missing_raw(self, adata_p3):
        adata_p3.raw = None
        issues = validate_p3(adata_p3)
        assert any("adata.raw" in i for i in issues)

    def test_alternative_reps(self, adata_p3):
        """X_pca should also be accepted."""
        del adata_p3.obsm["X_scvi"]
        adata_p3.obsm["X_pca"] = np.random.rand(adata_p3.n_obs, 50)
        assert validate_p3(adata_p3) == []

    def test_alternative_cluster_cols(self, adata_p3):
        """'cluster' column should also be accepted."""
        del adata_p3.obs["leiden"]
        adata_p3.obs["cluster"] = "0"
        assert validate_p3(adata_p3) == []


# ---------------------------------------------------------------------------
# validate_p35
# ---------------------------------------------------------------------------


class TestValidateP35:
    def test_valid(self, adata_p35):
        issues = validate_p35(adata_p35)
        assert issues == []

    def test_missing_signature_score(self, adata_p3):
        issues = validate_p35(adata_p3)
        assert any("signature_score" in i for i in issues)

    def test_missing_report(self, adata_p35):
        del adata_p35.uns["signature_score_report"]
        issues = validate_p35(adata_p35)
        assert any("signature_score_report" in i for i in issues)


# ---------------------------------------------------------------------------
# validate_p4
# ---------------------------------------------------------------------------


class TestValidateP4:
    def test_valid(self, adata_p4):
        issues = validate_p4(adata_p4)
        assert issues == []

    def test_missing_celltype(self, adata_p35):
        issues = validate_p4(adata_p35)
        assert any("celltype" in i and "celltype_broad" not in i for i in issues)

    def test_missing_celltype_broad(self, adata_p35):
        adata_p35.obs["celltype"] = "T cell"
        issues = validate_p4(adata_p35)
        assert any("celltype_broad" in i for i in issues)


# ---------------------------------------------------------------------------
# validate_checkpoint (dispatcher)
# ---------------------------------------------------------------------------


class TestValidateCheckpoint:
    def test_unknown_phase(self, adata_p1):
        with pytest.raises(ValueError, match="Unknown phase"):
            validate_checkpoint(adata_p1, phase="p99")

    def test_strict_raises(self, adata_p1):
        del adata_p1.obs["sample"]
        with pytest.raises(CheckpointValidationError):
            validate_checkpoint(adata_p1, phase="p1", strict=True)

    def test_non_strict_returns_issues(self, adata_p1):
        del adata_p1.obs["sample"]
        issues = validate_checkpoint(adata_p1, phase="p1", strict=False)
        assert len(issues) > 0

    def test_valid_returns_empty(self, adata_p1):
        issues = validate_checkpoint(adata_p1, phase="p1")
        assert issues == []

    # ── New slug-name tests ───────────────────────────────────────────────────

    def test_slug_qc_filter_accepted(self, adata_p1):
        """Slug 'qc_filter' is the preferred replacement for 'p1'."""
        issues = validate_checkpoint(adata_p1, phase="qc_filter")
        assert issues == []

    def test_slug_metadata_attach_accepted(self, adata_p2):
        issues = validate_checkpoint(adata_p2, phase="metadata_attach")
        assert issues == []

    def test_slug_preprocess_accepted(self, adata_p3):
        issues = validate_checkpoint(adata_p3, phase="preprocess")
        assert issues == []

    def test_slug_scoring_accepted(self, adata_p35):
        issues = validate_checkpoint(adata_p35, phase="scoring")
        assert issues == []

    def test_slug_celltype_manual_accepted(self, adata_p4):
        issues = validate_checkpoint(adata_p4, phase="celltype_manual")
        assert issues == []

    def test_legacy_p1_emits_deprecation_warning(self, adata_p1):
        """Old p1 code still works but emits DeprecationWarning."""
        import warnings

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            issues = validate_checkpoint(adata_p1, phase="p1")
            assert issues == []
            assert len(w) == 1
            assert issubclass(w[0].category, DeprecationWarning)
            assert "p1" in str(w[0].message)
            assert "qc_filter" in str(w[0].message)

    def test_legacy_p35_emits_deprecation_warning(self, adata_p35):
        """Old p35 code still works but emits DeprecationWarning."""
        import warnings

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            validate_checkpoint(adata_p35, phase="p35")
            assert len(w) == 1
            assert issubclass(w[0].category, DeprecationWarning)
            assert "scoring" in str(w[0].message)

    def test_unknown_slug_raises(self, adata_p1):
        """Unknown phase identifier (neither slug nor legacy code) raises ValueError."""
        with pytest.raises(ValueError, match="Unknown phase"):
            validate_checkpoint(adata_p1, phase="bogus_slug")


# ---------------------------------------------------------------------------
# validate_file
# ---------------------------------------------------------------------------


class TestValidateFile:
    def test_file_not_found_strict(self, tmp_path):
        with pytest.raises(CheckpointValidationError, match="File not found"):
            validate_file(tmp_path / "nonexistent.h5ad", phase="p1")

    def test_file_not_found_warn(self, tmp_path):
        issues = validate_file(tmp_path / "nonexistent.h5ad", phase="p1", strict=False)
        assert any("File not found" in i for i in issues)

    def test_roundtrip(self, adata_p1, tmp_path):
        path = tmp_path / "adata.raw.h5ad"
        adata_p1.write_h5ad(path)
        issues = validate_file(path, phase="p1")
        assert issues == []

    def test_fix_and_resave(self, adata_p1, tmp_path):
        del adata_p1.obs["raw_data_dir"]
        adata_p1.obs["batch"] = "batch_val"
        path = tmp_path / "adata.raw.h5ad"
        adata_p1.write_h5ad(path)

        issues = validate_file(path, phase="p1", fix=True)
        assert issues == []

        # Verify the fix persisted
        adata_reloaded = ad.read_h5ad(path)
        assert "raw_data_dir" in adata_reloaded.obs.columns
