"""Tests for sc_tools.pp.strategy -- Scale Strategy pattern."""

from __future__ import annotations

import numpy as np
import pytest
from anndata import AnnData


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
    adata = AnnData(X, obs={"library_id": (["L1"] * 100) + (["L2"] * 100)})
    adata.var_names = gene_names
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    return adata


class TestSmallStrategyPrepare:
    def test_prepare_sets_raw(self, adata_counts):
        from sc_tools.pp.strategy import SmallStrategy

        s = SmallStrategy()
        s.prepare(adata_counts)
        assert adata_counts.raw is not None

    def test_prepare_filters_genes(self, adata_counts):
        from sc_tools.pp.strategy import SmallStrategy

        s = SmallStrategy()
        n_before = adata_counts.n_vars
        s.prepare(adata_counts, filter_patterns=[r"^MT-"])
        assert adata_counts.n_vars == n_before - 3  # 3 MT genes removed


class TestSmallStrategySelectFeatures:
    def test_select_features_subsets_hvgs(self, adata_counts):
        import scanpy as sc

        from sc_tools.pp.strategy import SmallStrategy

        s = SmallStrategy()
        s.prepare(adata_counts)
        # Need to normalize first for HVG
        sc.pp.normalize_total(adata_counts)
        sc.pp.log1p(adata_counts)
        s.select_features(adata_counts, n_top_genes=100)
        assert adata_counts.n_vars == 100


class TestSmallStrategyReduceAndIntegrate:
    def test_none_integration(self, adata_counts):
        import scanpy as sc

        from sc_tools.pp.strategy import SmallStrategy

        s = SmallStrategy()
        s.prepare(adata_counts)
        sc.pp.normalize_total(adata_counts)
        sc.pp.log1p(adata_counts)
        s.select_features(adata_counts, n_top_genes=100)
        ctx = s.reduce_and_integrate(adata_counts, integration="none", n_comps=20)
        assert ctx is None
        assert "X_pca" in adata_counts.obsm
        assert adata_counts.obsm["X_pca"].shape[1] == 20

    def test_harmony_integration(self, adata_counts):
        try:
            import harmonypy  # noqa: F401
        except ImportError:
            pytest.skip("harmonypy not installed")
        import scanpy as sc

        from sc_tools.pp.strategy import SmallStrategy

        s = SmallStrategy()
        s.prepare(adata_counts)
        sc.pp.normalize_total(adata_counts)
        sc.pp.log1p(adata_counts)
        s.select_features(adata_counts, n_top_genes=100)
        ctx = s.reduce_and_integrate(
            adata_counts,
            integration="harmony",
            batch_key="library_id",
            n_comps=20,
        )
        assert ctx is None
        assert "X_pca_harmony" in adata_counts.obsm


class TestSmallStrategyEmbedAndCluster:
    def test_produces_leiden_and_umap(self, adata_counts):
        import scanpy as sc

        from sc_tools.pp.strategy import SmallStrategy

        s = SmallStrategy()
        s.prepare(adata_counts)
        sc.pp.normalize_total(adata_counts)
        sc.pp.log1p(adata_counts)
        s.select_features(adata_counts, n_top_genes=100)
        s.reduce_and_integrate(adata_counts, integration="none", n_comps=20)
        s.embed_and_cluster(adata_counts, resolution=0.5, n_neighbors=10)
        assert "leiden" in adata_counts.obs.columns
        assert "X_umap" in adata_counts.obsm


class TestSmallStrategyFullPipeline:
    def test_matches_visium_none_recipe(self, adata_counts):
        """Regression: SmallStrategy pipeline should produce same outputs as
        current _recipe_visium(integration='none')."""
        import scanpy as sc

        from sc_tools.pp.strategy import SmallStrategy

        s = SmallStrategy()
        s.prepare(adata_counts)
        sc.pp.normalize_total(adata_counts)
        sc.pp.log1p(adata_counts)
        s.select_features(adata_counts, n_top_genes=100)
        s.reduce_and_integrate(adata_counts, integration="none", n_comps=20)
        s.embed_and_cluster(adata_counts, resolution=0.5, n_neighbors=10)

        # All expected outputs present
        assert adata_counts.raw is not None
        assert "X_pca" in adata_counts.obsm
        assert "leiden" in adata_counts.obs.columns
        assert "X_umap" in adata_counts.obsm
        assert adata_counts.n_vars == 100


# ---------------------------------------------------------------------------
# LargeStrategy tests
# ---------------------------------------------------------------------------


class TestLargeStrategyInit:
    def test_name(self):
        from sc_tools.pp.strategy import LargeStrategy

        s = LargeStrategy()
        assert s.name == "large"

    def test_is_scale_strategy(self):
        from sc_tools.pp.strategy import LargeStrategy, ScaleStrategy

        s = LargeStrategy()
        assert isinstance(s, ScaleStrategy)

    def test_faiss_warning(self):
        """Should warn if FAISS unavailable (FAISS is not installed in test env)."""
        from sc_tools.pp.strategy import LargeStrategy

        with pytest.warns(UserWarning, match="FAISS"):
            LargeStrategy()


class TestLargeStrategyPrepare:
    def test_uses_layers_not_raw(self, adata_counts):
        from sc_tools.pp.strategy import LargeStrategy

        s = LargeStrategy()
        s.prepare(adata_counts)
        assert "raw_counts" in adata_counts.layers
        # Should NOT set .raw (memory saving)
        assert adata_counts.raw is None

    def test_filters_genes(self, adata_counts):
        from sc_tools.pp.strategy import LargeStrategy

        s = LargeStrategy()
        n_before = adata_counts.n_vars
        s.prepare(adata_counts, filter_patterns=[r"^MT-"])
        assert adata_counts.n_vars == n_before - 3


class TestLargeStrategySelectFeatures:
    def test_subsets_hvgs(self, adata_counts):
        import scanpy as sc

        from sc_tools.pp.strategy import LargeStrategy

        s = LargeStrategy()
        s.prepare(adata_counts)
        sc.pp.normalize_total(adata_counts)
        sc.pp.log1p(adata_counts)
        s.select_features(adata_counts, n_top_genes=100)
        assert adata_counts.n_vars == 100


class TestLargeStrategyReduceAndIntegrate:
    def test_none_returns_none(self, adata_counts):
        """integration='none': sparse PCA only, no subsampling."""
        import scanpy as sc

        from sc_tools.pp.strategy import LargeStrategy

        s = LargeStrategy()
        s.prepare(adata_counts)
        sc.pp.normalize_total(adata_counts)
        sc.pp.log1p(adata_counts)
        s.select_features(adata_counts, n_top_genes=100)
        ctx = s.reduce_and_integrate(adata_counts, integration="none", n_comps=20)
        assert ctx is None
        assert "X_pca" in adata_counts.obsm

    def test_harmony_returns_ctx(self, adata_counts):
        """Harmony: subsample -> correct -> project back."""
        try:
            import harmonypy  # noqa: F401
        except ImportError:
            pytest.skip("harmonypy not installed")
        import scanpy as sc

        from sc_tools.pp.strategy import LargeStrategy

        s = LargeStrategy(subsample_n=100)
        s.prepare(adata_counts)
        sc.pp.normalize_total(adata_counts)
        sc.pp.log1p(adata_counts)
        s.select_features(adata_counts, n_top_genes=100)
        ctx = s.reduce_and_integrate(
            adata_counts,
            integration="harmony",
            batch_key="library_id",
            n_comps=20,
        )
        assert ctx is not None
        assert ctx.subsample_idx is not None
        assert len(ctx.subsample_idx) == 100
        assert "X_pca_harmony" in adata_counts.obsm


class TestLargeStrategyEmbedAndCluster:
    def test_produces_outputs_with_projection(self, adata_counts):
        """Full pipeline: reduces, embeds, clusters, projects back."""
        import scanpy as sc

        from sc_tools.pp.strategy import LargeStrategy

        s = LargeStrategy(subsample_n=100, projection_k=10)
        s.prepare(adata_counts)
        sc.pp.normalize_total(adata_counts)
        sc.pp.log1p(adata_counts)
        s.select_features(adata_counts, n_top_genes=100)
        ctx = s.reduce_and_integrate(adata_counts, integration="none", n_comps=20)
        s.embed_and_cluster(adata_counts, ctx=ctx, resolution=0.5, n_neighbors=10)

        assert "leiden" in adata_counts.obs.columns
        assert "X_umap" in adata_counts.obsm
        assert "_sc_tools_is_representative" in adata_counts.obs.columns
        assert "_projection_info" in adata_counts.uns
        assert adata_counts.obs["_sc_tools_is_representative"].sum() == 100


# ---------------------------------------------------------------------------
# select_strategy() tests
# ---------------------------------------------------------------------------


class TestSelectStrategy:
    def test_small_adata_returns_small(self):
        """Small dataset should get SmallStrategy."""
        from sc_tools.pp.strategy import SmallStrategy, select_strategy

        adata = AnnData(np.random.randn(1000, 200).astype("float32"))
        strategy = select_strategy(adata)
        assert isinstance(strategy, SmallStrategy)

    def test_config_override_large(self):
        """Config backend='large' should force LargeStrategy."""
        from sc_tools.pp.strategy import LargeStrategy, select_strategy

        adata = AnnData(np.random.randn(100, 50).astype("float32"))
        strategy = select_strategy(adata, config={"backend": "large"})
        assert isinstance(strategy, LargeStrategy)

    def test_config_override_small(self):
        """Config backend='small' should force SmallStrategy."""
        from sc_tools.pp.strategy import SmallStrategy, select_strategy

        adata = AnnData(np.random.randn(100, 50).astype("float32"))
        strategy = select_strategy(adata, config={"backend": "small"})
        assert isinstance(strategy, SmallStrategy)

    def test_config_auto_defaults(self):
        """Config backend='auto' should auto-select (small for small data)."""
        from sc_tools.pp.strategy import SmallStrategy, select_strategy

        adata = AnnData(np.random.randn(100, 50).astype("float32"))
        strategy = select_strategy(adata, config={"backend": "auto"})
        assert isinstance(strategy, SmallStrategy)

    def test_mock_large_memory_triggers_large(self, monkeypatch):
        """When estimated peak exceeds 50% available memory -> LargeStrategy."""
        from sc_tools.pp import strategy as strat_mod
        from sc_tools.pp.strategy import LargeStrategy, select_strategy

        # Mock available memory to 4GB
        monkeypatch.setattr(strat_mod, "_get_available_memory_gb", lambda: 4.0)

        # Create adata where dense peak > 2GB (50% of 4GB)
        # n_obs * n_vars_est * 4 bytes / 1e9 > 2.0
        # 1M * 2000 * 4 = 8GB > 2GB
        adata = AnnData(np.random.randn(1_000_000, 10).astype("float32"))
        strategy = select_strategy(adata)
        assert isinstance(strategy, LargeStrategy)

    def test_platform_density_affects_threshold(self, monkeypatch):
        """IMC with 50 vars should be easier to keep as SmallStrategy."""
        from sc_tools.pp import strategy as strat_mod
        from sc_tools.pp.strategy import SmallStrategy, select_strategy

        monkeypatch.setattr(strat_mod, "_get_available_memory_gb", lambda: 8.0)
        # 500K cells * 50 vars * 4 bytes = 100MB -- well under threshold
        adata = AnnData(np.random.randn(500_000, 50).astype("float32"))
        strategy = select_strategy(adata, platform="imc")
        assert isinstance(strategy, SmallStrategy)

    def test_backed_adata_routes_to_large(self, tmp_path, monkeypatch):
        """Backed AnnData should always route to LargeStrategy."""
        from sc_tools.pp.strategy import LargeStrategy, select_strategy

        # Create a minimal backed AnnData
        h5_path = tmp_path / "test.h5ad"
        adata_mem = AnnData(np.random.randn(100, 50).astype("float32"))
        adata_mem.write_h5ad(h5_path)

        import anndata

        adata_backed = anndata.read_h5ad(h5_path, backed="r")
        strategy = select_strategy(adata_backed)
        assert isinstance(strategy, LargeStrategy)
        adata_backed.file.close()


class TestHasRapids:
    def test_returns_bool(self):
        from sc_tools.pp._gpu import has_rapids

        result = has_rapids()
        assert isinstance(result, bool)
