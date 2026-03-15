"""
Unit tests for sc_tools.qc.doublet: Solo-based doublet detection.

Tests use synthetic AnnData (50 cells, 100 genes, 2 batches).  scVI / SOLO
are mocked at the import site inside sc_tools.qc.doublet so the tests run
without a GPU and without real model training.  The scvi import-skip guard
is applied only to tests that actually call the mocked scvi symbols; the
modality-guard and import tests run unconditionally.
"""

from __future__ import annotations

import sys
import types
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
import scanpy as sc

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _make_sc_adata(n_obs: int = 50, n_vars: int = 100, n_batches: int = 2) -> sc.AnnData:
    """Synthetic single-cell AnnData with raw counts and two batches."""
    np.random.seed(0)
    X = np.random.negative_binomial(5, 0.3, (n_obs, n_vars)).astype(np.float32)
    obs = pd.DataFrame(
        {
            "library_id": [f"batch_{i % n_batches}" for i in range(n_obs)],
            "sample": [f"sample_{i % n_batches}" for i in range(n_obs)],
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])
    return sc.AnnData(X=X, obs=obs, var=var)


def _make_visium_adata(n_obs: int = 50, n_vars: int = 100) -> sc.AnnData:
    """Synthetic Visium (spot-based) AnnData."""
    adata = _make_sc_adata(n_obs=n_obs, n_vars=n_vars, n_batches=2)
    adata.obs["modality"] = "visium"
    return adata


def _make_imc_adata(n_obs: int = 50, n_vars: int = 40) -> sc.AnnData:
    """Synthetic IMC AnnData."""
    adata = _make_sc_adata(n_obs=n_obs, n_vars=n_vars, n_batches=2)
    adata.obs["modality"] = "imc"
    return adata


# ---------------------------------------------------------------------------
# Fake scvi stub injected into sys.modules before import
# ---------------------------------------------------------------------------


def _inject_scvi_stub():
    """
    Inject a minimal scvi stub into sys.modules so sc_tools.qc.doublet can be
    imported and its scvi symbols can be patched even when scvi is not installed.
    """
    if "scvi" in sys.modules:
        return  # real scvi already present

    scvi_mod = types.ModuleType("scvi")
    model_mod = types.ModuleType("scvi.model")
    external_mod = types.ModuleType("scvi.external")

    scvi_mod.model = model_mod
    scvi_mod.external = external_mod
    model_mod.SCVI = MagicMock(name="SCVI")
    external_mod.SOLO = MagicMock(name="SOLO")

    sys.modules["scvi"] = scvi_mod
    sys.modules["scvi.model"] = model_mod
    sys.modules["scvi.external"] = external_mod


_inject_scvi_stub()


# ---------------------------------------------------------------------------
# Helper: fake SOLO predict
# ---------------------------------------------------------------------------


def _mock_solo_predict(n_obs: int) -> pd.DataFrame:
    """Return a plausible SOLO predict DataFrame."""
    np.random.seed(1)
    scores = np.random.uniform(0, 1, n_obs)
    return pd.DataFrame({"doublet": scores, "singlet": 1 - scores})


def _build_solo_mock(n_obs: int) -> MagicMock:
    """Build a mock SOLO instance with a working predict method."""
    mock_instance = MagicMock()
    mock_instance.predict.return_value = _mock_solo_predict(n_obs)
    return mock_instance


# ---------------------------------------------------------------------------
# Import
# ---------------------------------------------------------------------------


def test_score_doublets_solo_importable():
    """score_doublets_solo must be importable from sc_tools.qc."""
    from sc_tools.qc import score_doublets_solo  # noqa: F401


# ---------------------------------------------------------------------------
# Output schema tests (use mocked scvi)
# ---------------------------------------------------------------------------


def test_score_doublets_solo_adds_obs_columns():
    """score_doublets_solo must add solo_doublet_score and solo_is_doublet to obs."""
    from sc_tools.qc import score_doublets_solo

    adata = _make_sc_adata(n_obs=50, n_vars=100)

    with (
        patch("sc_tools.qc.doublet.SCVI") as mock_scvi_cls,
        patch("sc_tools.qc.doublet.SOLO") as mock_solo_cls,
    ):
        mock_scvi_cls.return_value = MagicMock()
        mock_solo_cls.from_scvi_model.return_value = _build_solo_mock(adata.n_obs)

        result = score_doublets_solo(adata, batch_key="library_id", max_epochs=5, use_gpu=False)

    assert "solo_doublet_score" in result.obs.columns, "Missing solo_doublet_score in obs"
    assert "solo_is_doublet" in result.obs.columns, "Missing solo_is_doublet in obs"


def test_score_doublets_solo_score_dtype():
    """solo_doublet_score must be float; solo_is_doublet must be bool."""
    from sc_tools.qc import score_doublets_solo

    adata = _make_sc_adata(n_obs=50, n_vars=100)

    with (
        patch("sc_tools.qc.doublet.SCVI") as mock_scvi_cls,
        patch("sc_tools.qc.doublet.SOLO") as mock_solo_cls,
    ):
        mock_scvi_cls.return_value = MagicMock()
        mock_solo_cls.from_scvi_model.return_value = _build_solo_mock(adata.n_obs)

        result = score_doublets_solo(adata, batch_key="library_id", max_epochs=5, use_gpu=False)

    assert result.obs["solo_doublet_score"].dtype.kind == "f", "solo_doublet_score must be float"
    assert result.obs["solo_is_doublet"].dtype == bool, "solo_is_doublet must be bool"


def test_score_doublets_solo_uns_summary():
    """score_doublets_solo must populate uns['solo_summary'] with required keys."""
    from sc_tools.qc import score_doublets_solo

    adata = _make_sc_adata(n_obs=50, n_vars=100)

    with (
        patch("sc_tools.qc.doublet.SCVI") as mock_scvi_cls,
        patch("sc_tools.qc.doublet.SOLO") as mock_solo_cls,
    ):
        mock_scvi_cls.return_value = MagicMock()
        mock_solo_cls.from_scvi_model.return_value = _build_solo_mock(adata.n_obs)

        result = score_doublets_solo(adata, batch_key="library_id", max_epochs=5, use_gpu=False)

    assert "solo_summary" in result.uns, "Missing solo_summary in uns"
    summary = result.uns["solo_summary"]
    assert "n_doublets" in summary
    assert "pct_doublets" in summary
    assert "threshold" in summary
    assert isinstance(summary["n_doublets"], int)
    assert isinstance(summary["pct_doublets"], float)
    assert isinstance(summary["threshold"], float)


def test_score_doublets_solo_does_not_filter():
    """score_doublets_solo must NOT remove any cells -- only annotate."""
    from sc_tools.qc import score_doublets_solo

    adata = _make_sc_adata(n_obs=50, n_vars=100)
    n_obs = adata.n_obs

    with (
        patch("sc_tools.qc.doublet.SCVI") as mock_scvi_cls,
        patch("sc_tools.qc.doublet.SOLO") as mock_solo_cls,
    ):
        mock_scvi_cls.return_value = MagicMock()
        mock_solo_cls.from_scvi_model.return_value = _build_solo_mock(n_obs)

        result = score_doublets_solo(adata, batch_key="library_id", max_epochs=5, use_gpu=False)

    assert result.n_obs == n_obs, "score_doublets_solo must not filter any cells"


def test_score_doublets_solo_threshold_applied():
    """solo_is_doublet should be True exactly for scores above the threshold."""
    from sc_tools.qc import score_doublets_solo

    adata = _make_sc_adata(n_obs=50, n_vars=100)
    threshold = 0.5

    with (
        patch("sc_tools.qc.doublet.SCVI") as mock_scvi_cls,
        patch("sc_tools.qc.doublet.SOLO") as mock_solo_cls,
    ):
        mock_scvi_cls.return_value = MagicMock()
        mock_solo_cls.from_scvi_model.return_value = _build_solo_mock(adata.n_obs)

        result = score_doublets_solo(
            adata, batch_key="library_id", max_epochs=5, use_gpu=False, threshold=threshold
        )

    expected_doublets = result.obs["solo_doublet_score"] > threshold
    pd.testing.assert_series_equal(
        result.obs["solo_is_doublet"],
        expected_doublets,
        check_names=False,
    )


def test_score_doublets_solo_summary_counts():
    """n_doublets in summary must match the count of True values in solo_is_doublet."""
    from sc_tools.qc import score_doublets_solo

    adata = _make_sc_adata(n_obs=50, n_vars=100)

    with (
        patch("sc_tools.qc.doublet.SCVI") as mock_scvi_cls,
        patch("sc_tools.qc.doublet.SOLO") as mock_solo_cls,
    ):
        mock_scvi_cls.return_value = MagicMock()
        mock_solo_cls.from_scvi_model.return_value = _build_solo_mock(adata.n_obs)

        result = score_doublets_solo(adata, batch_key="library_id", max_epochs=5, use_gpu=False)

    expected_n = int(result.obs["solo_is_doublet"].sum())
    assert result.uns["solo_summary"]["n_doublets"] == expected_n


def test_score_doublets_solo_returns_anndata():
    """score_doublets_solo must return an AnnData object."""
    from sc_tools.qc import score_doublets_solo

    adata = _make_sc_adata(n_obs=50, n_vars=100)

    with (
        patch("sc_tools.qc.doublet.SCVI") as mock_scvi_cls,
        patch("sc_tools.qc.doublet.SOLO") as mock_solo_cls,
    ):
        mock_scvi_cls.return_value = MagicMock()
        mock_solo_cls.from_scvi_model.return_value = _build_solo_mock(adata.n_obs)

        result = score_doublets_solo(adata, batch_key="library_id", max_epochs=5, use_gpu=False)

    assert isinstance(result, sc.AnnData), "Return value must be AnnData"


# ---------------------------------------------------------------------------
# Modality guard tests (no scvi needed -- should raise before model training)
# ---------------------------------------------------------------------------


def test_score_doublets_solo_raises_for_visium_modality_obs():
    """ValueError must be raised when obs['modality'] == 'visium'."""
    from sc_tools.qc import score_doublets_solo

    adata = _make_visium_adata()
    with pytest.raises(ValueError, match="visium"):
        score_doublets_solo(adata, batch_key="library_id")


def test_score_doublets_solo_raises_for_imc_modality_obs():
    """ValueError must be raised when obs['modality'] == 'imc'."""
    from sc_tools.qc import score_doublets_solo

    adata = _make_imc_adata()
    with pytest.raises(ValueError, match="imc"):
        score_doublets_solo(adata, batch_key="library_id")


def test_score_doublets_solo_raises_via_modality_param():
    """ValueError must be raised when modality parameter is 'visium'."""
    from sc_tools.qc import score_doublets_solo

    adata = _make_sc_adata()
    with pytest.raises(ValueError, match="visium"):
        score_doublets_solo(adata, batch_key="library_id", modality="visium")


def test_score_doublets_solo_raises_via_modality_param_imc():
    """ValueError must be raised when modality parameter is 'imc'."""
    from sc_tools.qc import score_doublets_solo

    adata = _make_sc_adata()
    with pytest.raises(ValueError, match="imc"):
        score_doublets_solo(adata, batch_key="library_id", modality="imc")
