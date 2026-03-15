"""Tests for run_resolvi() integration method.

Follows the TDD pattern from test_pp.py: skip if scvi-tools is not installed,
verify ImportError path when missing, and run a minimal training check when available.
"""

from __future__ import annotations

import numpy as np
import pytest
from anndata import AnnData
from scipy import sparse

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _resolvi_available() -> bool:
    try:
        from scvi.external import RESOLVI  # noqa: F401

        return True
    except ImportError:
        return False


def _make_small_spatial_adata(
    n_obs: int = 100,
    n_vars: int = 50,
    n_batches: int = 2,
    seed: int = 0,
) -> AnnData:
    """Synthetic AnnData with raw counts, spatial coords, and batch labels."""
    rng = np.random.default_rng(seed)
    X = rng.negative_binomial(5, 0.3, (n_obs, n_vars)).astype("float32")

    batch_labels = [f"batch_{i % n_batches}" for i in range(n_obs)]
    spatial = rng.uniform(0, 1000, (n_obs, 2)).astype("float32")

    adata = AnnData(
        X=sparse.csr_matrix(X),
        obs={"library_id": batch_labels},
    )
    adata.var_names = [f"GENE{i}" for i in range(n_vars)]
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adata.obsm["spatial"] = spatial
    # RESOLVI requires obsm['X_spatial']
    adata.obsm["X_spatial"] = spatial.copy()
    return adata


# ---------------------------------------------------------------------------
# Import-error path (when scvi-tools is absent)
# ---------------------------------------------------------------------------


class TestRunResolviImportError:
    def test_raises_import_error_when_scvi_missing(self):
        """run_resolvi should raise ImportError with install hint if scvi-tools absent."""
        if _resolvi_available():
            pytest.skip("scvi-tools is installed; cannot test ImportError path")

        from sc_tools.pp import run_resolvi

        adata = _make_small_spatial_adata()
        with pytest.raises(ImportError, match="scvi-tools"):
            run_resolvi(adata)


# ---------------------------------------------------------------------------
# Functional tests (skipped when scvi-tools absent)
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not _resolvi_available(), reason="scvi-tools / RESOLVI not installed")
class TestRunResolvi:
    def test_stores_latent_in_obsm(self):
        """run_resolvi must write X_resolvi into adata.obsm."""
        from sc_tools.pp import run_resolvi

        adata = _make_small_spatial_adata()
        run_resolvi(
            adata, batch_key="library_id", max_epochs=2, use_gpu=False, downsample_counts=False
        )
        assert "X_resolvi" in adata.obsm, "X_resolvi not written to adata.obsm"

    def test_latent_shape_matches_n_obs(self):
        """Latent embedding must have shape (n_obs, n_latent)."""
        from sc_tools.pp import run_resolvi

        n_obs, n_latent = 100, 5
        adata = _make_small_spatial_adata(n_obs=n_obs)
        run_resolvi(
            adata,
            batch_key="library_id",
            n_latent=n_latent,
            max_epochs=2,
            use_gpu=False,
            downsample_counts=False,
        )
        emb = adata.obsm["X_resolvi"]
        assert emb.shape == (n_obs, n_latent), (
            f"Expected shape ({n_obs}, {n_latent}), got {emb.shape}"
        )

    def test_returns_trained_model(self):
        """run_resolvi must return the trained model (not None)."""
        from sc_tools.pp import run_resolvi

        adata = _make_small_spatial_adata()
        model = run_resolvi(
            adata, batch_key="library_id", max_epochs=2, use_gpu=False, downsample_counts=False
        )
        assert model is not None, "run_resolvi returned None; expected a trained model"

    def test_layer_argument_respected(self):
        """run_resolvi must use the specified layer when provided."""
        from sc_tools.pp import run_resolvi

        adata = _make_small_spatial_adata()
        # Normalize X but put raw counts in a layer
        adata.layers["counts"] = adata.X.copy()
        import scanpy as sc

        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        model = run_resolvi(
            adata,
            batch_key="library_id",
            layer="counts",
            max_epochs=2,
            use_gpu=False,
            downsample_counts=False,
        )
        assert "X_resolvi" in adata.obsm
        assert model is not None

    def test_spatial_key_auto_copied(self):
        """run_resolvi must copy obsm['spatial'] -> obsm['X_spatial'] when absent."""
        from sc_tools.pp import run_resolvi

        adata = _make_small_spatial_adata()
        # Remove X_spatial so the function must create it
        del adata.obsm["X_spatial"]
        assert "X_spatial" not in adata.obsm

        run_resolvi(
            adata, batch_key="library_id", max_epochs=2, use_gpu=False, downsample_counts=False
        )
        assert "X_resolvi" in adata.obsm

    def test_params_stored_in_uns(self):
        """run_resolvi should record params in adata.uns['resolvi_params']."""
        from sc_tools.pp import run_resolvi

        adata = _make_small_spatial_adata()
        run_resolvi(
            adata, batch_key="library_id", max_epochs=2, use_gpu=False, downsample_counts=False
        )
        assert "resolvi_params" in adata.uns
        params = adata.uns["resolvi_params"]
        assert params["batch_key"] == "library_id"
        assert "max_epochs" in params

    def test_exported_from_pp(self):
        """run_resolvi must be importable directly from sc_tools.pp."""
        from sc_tools import pp

        assert hasattr(pp, "run_resolvi"), "run_resolvi not exported from sc_tools.pp"

    def test_in_transcriptomic_methods(self):
        """resolvi must appear in the transcriptomic integration candidate list."""
        from sc_tools.bm.integration import _TRANSCRIPTOMIC_METHODS

        assert "resolvi" in _TRANSCRIPTOMIC_METHODS, (
            "'resolvi' not in _TRANSCRIPTOMIC_METHODS benchmark candidates"
        )
