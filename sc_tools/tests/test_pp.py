"""Unit tests for sc_tools.pp preprocessing module.

Tests use synthetic AnnData fixtures. Integration tests (scVI, Harmony, CytoVI)
are skipped if the required packages are not installed.
"""

from __future__ import annotations

import numpy as np
import pytest
from anndata import AnnData
from scipy import sparse

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def adata_counts():
    """Synthetic count matrix with MT/RP/HB genes and library_id."""
    np.random.seed(42)
    n_obs, n_vars = 200, 500
    X = np.random.negative_binomial(5, 0.3, (n_obs, n_vars)).astype("float32")

    gene_names = [f"GENE{i}" for i in range(n_vars - 8)]
    gene_names += ["MT-CO1", "MT-CO2", "MT-ND1"]
    gene_names += ["RPS6", "RPL11"]
    gene_names += ["HBA1", "HBA2", "HBB"]

    adata = AnnData(
        X,
        obs={"library_id": (["L1"] * 100) + (["L2"] * 100)},
    )
    adata.var_names = gene_names
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    return adata


@pytest.fixture()
def adata_imc():
    """Synthetic IMC (protein) data with 40 markers."""
    np.random.seed(123)
    n_obs, n_vars = 150, 40
    X = np.random.poisson(10, (n_obs, n_vars)).astype("float32")
    adata = AnnData(
        X,
        obs={"library_id": (["ROI1"] * 75) + (["ROI2"] * 75)},
    )
    adata.var_names = [f"Marker{i}" for i in range(n_vars)]
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    return adata


@pytest.fixture()
def adata_sparse(adata_counts):
    """Same as adata_counts but with sparse X."""
    adata = adata_counts.copy()
    adata.X = sparse.csr_matrix(adata.X)
    return adata


# ---------------------------------------------------------------------------
# GPU backend
# ---------------------------------------------------------------------------


class TestGPUBackend:
    def test_get_backend_returns_scanpy(self):
        from sc_tools.pp._gpu import get_backend

        backend, name = get_backend()
        # On most test environments, rapids won't be installed
        assert name in ("scanpy", "rapids")
        assert hasattr(backend, "pp")


# ---------------------------------------------------------------------------
# Normalization
# ---------------------------------------------------------------------------


class TestBackupRaw:
    def test_backup_creates_raw(self, adata_counts):
        from sc_tools.pp import backup_raw

        assert adata_counts.raw is None
        backup_raw(adata_counts)
        assert adata_counts.raw is not None
        assert adata_counts.raw.n_vars == 500

    def test_backup_noop_if_exists(self, adata_counts):
        from sc_tools.pp import backup_raw

        adata_counts.raw = adata_counts.copy()
        original_raw_shape = adata_counts.raw.shape
        backup_raw(adata_counts)
        assert adata_counts.raw.shape == original_raw_shape


class TestNormalizeTotal:
    def test_normalize_inplace(self, adata_counts):
        from sc_tools.pp import normalize_total

        original_sums = adata_counts.X.sum(axis=1)
        normalize_total(adata_counts, target_sum=1e4)
        new_sums = adata_counts.X.sum(axis=1)
        np.testing.assert_allclose(new_sums, 1e4, rtol=1e-5)
        assert not np.allclose(original_sums, new_sums)

    def test_normalize_sparse(self, adata_sparse):
        from sc_tools.pp import normalize_total

        normalize_total(adata_sparse, target_sum=1e4)
        sums = np.asarray(adata_sparse.X.sum(axis=1)).flatten()
        np.testing.assert_allclose(sums, 1e4, rtol=1e-5)


class TestLogTransform:
    def test_log_transform(self, adata_counts):
        from sc_tools.pp import log_transform, normalize_total

        normalize_total(adata_counts)
        before = adata_counts.X.copy()
        log_transform(adata_counts)
        # log1p should reduce values
        assert adata_counts.X.max() < before.max()
        # Check log1p identity
        np.testing.assert_allclose(adata_counts.X, np.log1p(before), rtol=1e-5)


class TestScale:
    def test_scale_clamps(self, adata_counts):
        from sc_tools.pp import log_transform, normalize_total, scale

        normalize_total(adata_counts)
        log_transform(adata_counts)
        scale(adata_counts, max_value=10)
        assert adata_counts.X.max() <= 10.0 + 1e-6


class TestArcsinhTransform:
    def test_arcsinh_inplace(self, adata_imc):
        from sc_tools.pp import arcsinh_transform

        original = adata_imc.X.copy()
        arcsinh_transform(adata_imc, cofactor=5)
        expected = np.arcsinh(original / 5)
        np.testing.assert_allclose(adata_imc.X, expected, rtol=1e-5)

    def test_arcsinh_copy(self, adata_imc):
        from sc_tools.pp import arcsinh_transform

        original = adata_imc.X.copy()
        result = arcsinh_transform(adata_imc, cofactor=5, inplace=False)
        assert result is not None
        # Original unchanged
        np.testing.assert_array_equal(adata_imc.X, original)
        # Result is transformed
        np.testing.assert_allclose(result.X, np.arcsinh(original / 5), rtol=1e-5)

    def test_arcsinh_sparse(self, adata_imc):
        from sc_tools.pp import arcsinh_transform

        adata_imc.X = sparse.csr_matrix(adata_imc.X)
        arcsinh_transform(adata_imc, cofactor=5)
        # Should be dense after arcsinh
        assert not sparse.issparse(adata_imc.X)


class TestFilterGenesByPattern:
    def test_removes_mt_rp_hb(self, adata_counts):
        from sc_tools.pp import filter_genes_by_pattern

        n_before = adata_counts.n_vars
        filter_genes_by_pattern(adata_counts)
        n_after = adata_counts.n_vars
        assert n_after == n_before - 8  # 3 MT + 2 RP + 3 HB
        # Verify none of the pattern genes remain
        remaining = set(adata_counts.var_names)
        for g in ["MT-CO1", "MT-CO2", "MT-ND1", "RPS6", "RPL11", "HBA1", "HBA2", "HBB"]:
            assert g not in remaining

    def test_custom_patterns(self, adata_counts):
        from sc_tools.pp import filter_genes_by_pattern

        filter_genes_by_pattern(adata_counts, patterns=[r"^MT-"])
        remaining = set(adata_counts.var_names)
        assert "MT-CO1" not in remaining
        assert "RPS6" in remaining  # RP genes not removed

    def test_keep_mode(self, adata_counts):
        from sc_tools.pp import filter_genes_by_pattern

        filter_genes_by_pattern(adata_counts, patterns=[r"^MT-"], exclude=False)
        assert adata_counts.n_vars == 3
        assert all(g.startswith("MT-") for g in adata_counts.var_names)

    def test_case_insensitive(self):
        from sc_tools.pp import filter_genes_by_pattern

        adata = AnnData(np.ones((5, 4), dtype="float32"))
        adata.var_names = ["mt-co1", "Mt-Nd1", "GENE1", "GENE2"]
        filter_genes_by_pattern(adata, patterns=[r"^MT-"])
        assert adata.n_vars == 2


# ---------------------------------------------------------------------------
# Dimensionality reduction & clustering
# ---------------------------------------------------------------------------


class TestPCA:
    def test_pca_basic(self, adata_counts):
        from sc_tools.pp import log_transform, normalize_total, pca

        normalize_total(adata_counts)
        log_transform(adata_counts)
        pca(adata_counts, n_comps=20)
        assert "X_pca" in adata_counts.obsm
        assert adata_counts.obsm["X_pca"].shape == (200, 20)

    def test_pca_clamps_n_comps(self):
        from sc_tools.pp import pca

        # Small adata where n_comps > min(n_obs, n_vars) - 1
        adata = AnnData(np.random.randn(10, 5).astype("float32"))
        pca(adata, n_comps=50)
        assert "X_pca" in adata.obsm
        assert adata.obsm["X_pca"].shape[1] <= 4  # min(10, 5) - 1


class TestNeighbors:
    def test_auto_detect_use_rep(self, adata_counts):
        from sc_tools.pp import log_transform, neighbors, normalize_total, pca

        normalize_total(adata_counts)
        log_transform(adata_counts)
        pca(adata_counts, n_comps=20)
        neighbors(adata_counts, n_neighbors=10)
        assert "neighbors" in adata_counts.uns

    def test_explicit_use_rep(self, adata_counts):
        from sc_tools.pp import neighbors

        # Fake a scVI latent
        adata_counts.obsm["X_scVI"] = np.random.randn(200, 10).astype("float32")
        neighbors(adata_counts, use_rep="X_scVI", n_neighbors=10)
        assert "neighbors" in adata_counts.uns


class TestLeiden:
    def test_leiden_basic(self, adata_counts):
        from sc_tools.pp import leiden, log_transform, neighbors, normalize_total, pca

        normalize_total(adata_counts)
        log_transform(adata_counts)
        pca(adata_counts, n_comps=20)
        neighbors(adata_counts, n_neighbors=10)
        leiden(adata_counts, resolution=0.5)
        assert "leiden" in adata_counts.obs.columns

    def test_leiden_custom_key(self, adata_counts):
        from sc_tools.pp import leiden, log_transform, neighbors, normalize_total, pca

        normalize_total(adata_counts)
        log_transform(adata_counts)
        pca(adata_counts, n_comps=20)
        neighbors(adata_counts, n_neighbors=10)
        leiden(adata_counts, resolution=0.5, key_added="my_clusters")
        assert "my_clusters" in adata_counts.obs.columns


class TestCluster:
    def test_cluster_convenience(self, adata_counts):
        from sc_tools.pp import cluster, log_transform, normalize_total, pca

        normalize_total(adata_counts)
        log_transform(adata_counts)
        pca(adata_counts, n_comps=20)
        cluster(adata_counts, resolution=0.5, n_neighbors=10)
        assert "leiden" in adata_counts.obs.columns
        assert "X_umap" in adata_counts.obsm

    def test_cluster_no_umap(self, adata_counts):
        from sc_tools.pp import cluster, log_transform, normalize_total, pca

        normalize_total(adata_counts)
        log_transform(adata_counts)
        pca(adata_counts, n_comps=20)
        cluster(adata_counts, resolution=0.5, n_neighbors=10, run_umap=False)
        assert "leiden" in adata_counts.obs.columns
        assert "X_umap" not in adata_counts.obsm


class TestUMAP:
    def test_umap(self, adata_counts):
        from sc_tools.pp import log_transform, neighbors, normalize_total, pca, umap

        normalize_total(adata_counts)
        log_transform(adata_counts)
        pca(adata_counts, n_comps=20)
        neighbors(adata_counts, n_neighbors=10)
        umap(adata_counts)
        assert "X_umap" in adata_counts.obsm
        assert adata_counts.obsm["X_umap"].shape == (200, 2)


# ---------------------------------------------------------------------------
# UTAG (soft dependency)
# ---------------------------------------------------------------------------


class TestRunUTAG:
    def test_utag_import_error(self, adata_counts):
        """Should raise ImportError with install hint if utag not installed."""
        from sc_tools.pp import run_utag

        # UTAG is unlikely to be installed in test env
        try:
            import utag  # noqa: F401

            pytest.skip("utag is installed; cannot test ImportError")
        except ImportError:
            with pytest.raises(ImportError, match="UTAG is required"):
                run_utag(adata_counts)


# ---------------------------------------------------------------------------
# Integration (soft dependencies)
# ---------------------------------------------------------------------------


def _scvi_available() -> bool:
    try:
        import scvi  # noqa: F401

        return True
    except ImportError:
        return False


class TestRunScVI:
    def test_scvi_import_error(self):
        """Should raise ImportError with install hint if scvi not installed."""
        try:
            import scvi  # noqa: F401

            pytest.skip("scvi-tools is installed; skip import error test")
        except ImportError:
            from sc_tools.pp import run_scvi

            adata = AnnData(np.ones((10, 5), dtype="float32"))
            with pytest.raises(ImportError, match="scvi-tools"):
                run_scvi(adata)

    @pytest.mark.skipif(
        not _scvi_available(),
        reason="scvi-tools not installed",
    )
    def test_scvi_runs(self, adata_counts):
        from sc_tools.pp import run_scvi

        run_scvi(adata_counts, batch_key="library_id", max_epochs=5, n_latent=5)
        assert "X_scVI" in adata_counts.obsm
        assert "scvi_params" in adata_counts.uns


class TestRunHarmony:
    def test_harmony_import_error(self):
        """Should raise ImportError with install hint if harmonypy not installed."""
        try:
            import harmonypy  # noqa: F401

            pytest.skip("harmonypy is installed; skip import error test")
        except ImportError:
            from sc_tools.pp import run_harmony

            adata = AnnData(np.ones((10, 5), dtype="float32"))
            with pytest.raises(ImportError, match="harmonypy"):
                run_harmony(adata)

    def test_harmony_missing_basis(self, adata_counts):
        from sc_tools.pp import run_harmony

        try:
            import harmonypy  # noqa: F401
        except ImportError:
            pytest.skip("harmonypy not installed")

        with pytest.raises(ValueError, match="not found in adata.obsm"):
            run_harmony(adata_counts, basis="X_pca")


class TestRunCytoVI:
    def test_cytovi_import_error(self):
        """Should raise ImportError with install hint if scvi not installed."""
        try:
            import scvi  # noqa: F401

            pytest.skip("scvi-tools is installed; skip import error test")
        except ImportError:
            from sc_tools.pp import run_cytovi

            adata = AnnData(np.ones((10, 5), dtype="float32"))
            with pytest.raises(ImportError, match="scvi-tools"):
                run_cytovi(adata)


# ---------------------------------------------------------------------------
# Recipes
# ---------------------------------------------------------------------------


class TestPreprocessRecipe:
    def test_invalid_modality(self, adata_counts):
        from sc_tools.pp import preprocess

        with pytest.raises(ValueError, match="Unknown modality"):
            preprocess(adata_counts, modality="invalid")

    def test_invalid_integration(self, adata_counts):
        from sc_tools.pp import preprocess

        with pytest.raises(ValueError, match="Unknown integration"):
            preprocess(adata_counts, integration="invalid")

    def test_visium_no_integration(self, adata_counts):
        from sc_tools.pp import preprocess

        result = preprocess(
            adata_counts,
            modality="visium",
            integration="none",
            n_top_genes=100,
            resolution=0.5,
            copy=True,
        )
        assert result.raw is not None
        assert "leiden" in result.obs.columns
        assert "X_umap" in result.obsm
        assert "X_pca" in result.obsm
        # MT/RP/HB genes should be filtered
        assert result.n_vars <= 100

    def test_xenium_recipe(self, adata_counts):
        from sc_tools.pp import preprocess

        result = preprocess(
            adata_counts,
            modality="xenium",
            integration="none",
            n_top_genes=100,
            resolution=0.5,
            copy=True,
        )
        assert result.raw is not None
        assert "leiden" in result.obs.columns
        assert "X_umap" in result.obsm

    def test_cosmx_recipe(self, adata_counts):
        from sc_tools.pp import preprocess

        result = preprocess(
            adata_counts,
            modality="cosmx",
            integration="none",
            n_top_genes=100,
            resolution=0.5,
            copy=True,
        )
        assert result.raw is not None
        assert "leiden" in result.obs.columns

    def test_imc_recipe(self, adata_imc):
        from sc_tools.pp import preprocess

        result = preprocess(
            adata_imc,
            modality="imc",
            integration="none",
            resolution=0.5,
            copy=True,
        )
        assert result.raw is not None
        assert "leiden" in result.obs.columns
        assert "X_pca" in result.obsm
        # After arcsinh + scale, PCA components should be present
        assert result.obsm["X_pca"].shape[1] <= result.n_vars

    def test_copy_preserves_original(self, adata_counts):
        from sc_tools.pp import preprocess

        original_shape = adata_counts.shape
        _ = preprocess(
            adata_counts,
            modality="visium",
            integration="none",
            n_top_genes=100,
            copy=True,
        )
        # Original should be unchanged
        assert adata_counts.shape == original_shape
        assert adata_counts.raw is None

    def test_inplace_modifies(self, adata_counts):
        from sc_tools.pp import preprocess

        preprocess(
            adata_counts,
            modality="visium",
            integration="none",
            n_top_genes=100,
        )
        # Original should be modified
        assert adata_counts.raw is not None
        assert "leiden" in adata_counts.obs.columns


class TestRecipeUsesStrategy:
    """Verify recipes accept and use strategy parameter."""

    def test_preprocess_uses_strategy_internally(self, adata_counts):
        """preprocess() should work exactly as before (behavioral regression test)."""
        from sc_tools.pp import preprocess

        result = preprocess(
            adata_counts,
            modality="visium",
            integration="none",
            n_top_genes=100,
            resolution=0.5,
            copy=True,
        )
        assert result.raw is not None
        assert "leiden" in result.obs.columns
        assert "X_umap" in result.obsm
        assert "X_pca" in result.obsm

    def test_recipe_accepts_strategy_kwarg(self, adata_counts):
        """Recipes should accept an explicit strategy parameter."""
        from sc_tools.pp import preprocess
        from sc_tools.pp.strategy import SmallStrategy

        result = preprocess(
            adata_counts,
            modality="visium",
            integration="none",
            n_top_genes=100,
            resolution=0.5,
            strategy=SmallStrategy(),
            copy=True,
        )
        assert result.raw is not None
        assert "leiden" in result.obs.columns

    def test_targeted_panel_xenium(self, adata_counts):
        """Xenium should still work after merge into targeted panel."""
        from sc_tools.pp import preprocess

        result = preprocess(
            adata_counts,
            modality="xenium",
            integration="none",
            n_top_genes=100,
            resolution=0.5,
            copy=True,
        )
        assert result.raw is not None
        assert "leiden" in result.obs.columns

    def test_targeted_panel_cosmx(self, adata_counts):
        """CosMx should still work after merge into targeted panel."""
        from sc_tools.pp import preprocess

        result = preprocess(
            adata_counts,
            modality="cosmx",
            integration="none",
            n_top_genes=100,
            resolution=0.5,
            copy=True,
        )
        assert result.raw is not None
        assert "leiden" in result.obs.columns

    def test_imc_recipe_with_strategy(self, adata_imc):
        """IMC recipe should work with strategy pattern."""
        from sc_tools.pp import preprocess

        result = preprocess(
            adata_imc,
            modality="imc",
            integration="none",
            resolution=0.5,
            copy=True,
        )
        assert result.raw is not None
        assert "leiden" in result.obs.columns
        assert "X_pca" in result.obsm


class TestAutoUseRep:
    def test_priority_order(self):
        from sc_tools.pp.reduce import _auto_use_rep

        adata = AnnData(np.ones((5, 3), dtype="float32"))

        # No embeddings -> default X_pca
        assert _auto_use_rep(adata, None) == "X_pca"

        # Add X_pca
        adata.obsm["X_pca"] = np.ones((5, 2))
        assert _auto_use_rep(adata, None) == "X_pca"

        # Add X_pca_harmony (higher priority)
        adata.obsm["X_pca_harmony"] = np.ones((5, 2))
        assert _auto_use_rep(adata, None) == "X_pca_harmony"

        # Add X_scVI (highest priority)
        adata.obsm["X_scVI"] = np.ones((5, 2))
        assert _auto_use_rep(adata, None) == "X_scVI"

        # Explicit override
        assert _auto_use_rep(adata, "X_custom") == "X_custom"


# ---------------------------------------------------------------------------
# TestTargetedPanelScVI (BM-07)
# ---------------------------------------------------------------------------


class TestTargetedPanelScVI:
    def test_scvi_preserves_raw_counts(self):
        """_recipe_targeted_panel with scVI should NOT normalize X."""
        from unittest.mock import MagicMock, patch

        from sc_tools.pp.recipes import _recipe_targeted_panel
        from sc_tools.pp.strategy import SmallStrategy, SubsampleContext

        rng = np.random.RandomState(42)
        n_obs, n_vars = 100, 50
        X = rng.negative_binomial(5, 0.3, (n_obs, n_vars)).astype(np.float32)
        adata = AnnData(X=X.copy())
        adata.obs["batch"] = [f"batch_{i % 3}" for i in range(n_obs)]

        strategy = SmallStrategy()
        ctx = SubsampleContext()

        with patch.object(strategy, "reduce_and_integrate", return_value=ctx), \
             patch.object(strategy, "embed_and_cluster"):
            _recipe_targeted_panel(
                adata,
                batch_key="batch",
                integration="scvi",
                n_top_genes=20,
                resolution=0.5,
                filter_patterns=None,
                use_gpu=False,
                strategy=strategy,
            )

        # Raw counts should be preserved (non-negative integers)
        assert np.all(adata.X >= 0), "X has negative values after scVI path"
        # Check values are close to integers (raw counts are integers)
        assert np.allclose(adata.X[adata.X > 0], np.round(adata.X[adata.X > 0]), atol=0.01), (
            "X values are not integers after scVI path -- normalization ran incorrectly"
        )

    def test_harmony_normalizes(self):
        """_recipe_targeted_panel with harmony should normalize X."""
        from unittest.mock import patch

        from sc_tools.pp.recipes import _recipe_targeted_panel
        from sc_tools.pp.strategy import SmallStrategy, SubsampleContext

        rng = np.random.RandomState(42)
        n_obs, n_vars = 100, 50
        X = rng.negative_binomial(5, 0.3, (n_obs, n_vars)).astype(np.float32)
        adata = AnnData(X=X.copy())
        adata.obs["batch"] = [f"batch_{i % 3}" for i in range(n_obs)]
        original_X = adata.X.copy()

        strategy = SmallStrategy()
        ctx = SubsampleContext()

        with patch.object(strategy, "reduce_and_integrate", return_value=ctx), \
             patch.object(strategy, "embed_and_cluster"):
            _recipe_targeted_panel(
                adata,
                batch_key="batch",
                integration="harmony",
                n_top_genes=20,
                resolution=0.5,
                filter_patterns=None,
                use_gpu=False,
                strategy=strategy,
            )

        # X should have been normalized (log1p -> non-integer values)
        assert not np.array_equal(adata.X, original_X), "X unchanged after harmony path"

    def test_scvi_uses_seurat_v3_features(self):
        """scVI path should call select_features with flavor='seurat_v3'."""
        from unittest.mock import MagicMock, patch

        from sc_tools.pp.recipes import _recipe_targeted_panel
        from sc_tools.pp.strategy import SmallStrategy, SubsampleContext

        rng = np.random.RandomState(42)
        n_obs, n_vars = 100, 50
        X = rng.negative_binomial(5, 0.3, (n_obs, n_vars)).astype(np.float32)
        adata = AnnData(X=X.copy())
        adata.obs["batch"] = [f"batch_{i % 3}" for i in range(n_obs)]

        strategy = SmallStrategy()
        ctx = SubsampleContext()

        with patch.object(strategy, "select_features") as mock_select, \
             patch.object(strategy, "reduce_and_integrate", return_value=ctx), \
             patch.object(strategy, "embed_and_cluster"):
            _recipe_targeted_panel(
                adata,
                batch_key="batch",
                integration="scvi",
                n_top_genes=20,
                resolution=0.5,
                filter_patterns=None,
                use_gpu=False,
                strategy=strategy,
            )

        mock_select.assert_called_once()
        call_kwargs = mock_select.call_args
        assert call_kwargs[1].get("flavor") == "seurat_v3" or \
               (len(call_kwargs[0]) > 2 and call_kwargs[0][2] == "seurat_v3"), \
            f"select_features not called with flavor='seurat_v3': {call_kwargs}"
