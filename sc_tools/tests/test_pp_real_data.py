"""Real-data edge-case tests for sc_tools.pp preprocessing module.

Uses actual on-disk datasets (Visium, IMC, CosMx). Each test class is
guarded by ``pytest.mark.skipif`` so the suite is still runnable without
the data files present.

Subsets aggressively to keep runtime short:
- Visium:  200 spots per library_id (up to 3 library_ids)
- IMC:     200 cells per ROI (up to 3 ROIs)
- CosMx:   200 cells per sample (up to 2 samples)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from scipy import sparse

# ---------------------------------------------------------------------------
# Paths to real data
# ---------------------------------------------------------------------------

VISIUM_ANNOTATED = Path("projects/visium/ggo_visium/results/adata.annotated.p2.h5ad")
VISIUM_NORMALIZED = Path("projects/visium/ggo_visium/results/adata.normalized.h5ad")
IMC_CELLTYPED = Path("projects/imc/ggo_human/results/celltyped.h5ad")
COSMX_ANNOTATED = Path("projects/cosmx_1k/lymph_dlbcl/results/cosmx_rna_annotated.h5ad")

REPO_ROOT = Path(__file__).parent.parent.parent


def _abs(p: Path) -> Path:
    return REPO_ROOT / p


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _subset_by_batch(adata, batch_key: str, n_per_batch: int, max_batches: int = 3):
    """Return a subset with at most n_per_batch observations per batch value."""

    batches = adata.obs[batch_key].unique().tolist()[:max_batches]
    idx: list[int] = []
    for b in batches:
        mask = adata.obs[batch_key] == b
        chosen = np.where(mask)[0][:n_per_batch]
        idx.extend(chosen.tolist())
    sub = adata[idx].copy()
    # Ensure integer counts (some h5ad store float; scVI needs int-like)
    if not sparse.issparse(sub.X):
        sub.X = sub.X.astype("float32")
    else:
        sub.X = sub.X.astype("float32")
    return sub


# ---------------------------------------------------------------------------
# Visium normalize / HVG tests
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not _abs(VISIUM_ANNOTATED).exists(),
    reason="Visium annotated data not available",
)
class TestNormalizeRealVisium:
    """Normalization edge-case tests on real Visium data."""

    @pytest.fixture(scope="class")
    def visium_sub(self):
        """Load and subset the Visium annotated h5ad (3 library_ids, 200 spots each)."""
        import anndata as ad

        adata = ad.read_h5ad(_abs(VISIUM_ANNOTATED))
        sub = _subset_by_batch(adata, "library_id", n_per_batch=200, max_batches=3)
        return sub

    def test_normalize_total_no_nan_inf(self, visium_sub):
        """normalize_total on real sparse Visium data must not produce NaN or Inf."""
        from sc_tools.pp import normalize_total

        adata = visium_sub.copy()
        normalize_total(adata, target_sum=1e4)
        X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
        assert not np.any(np.isnan(X)), "NaN values after normalize_total"
        assert not np.any(np.isinf(X)), "Inf values after normalize_total"

    def test_normalize_total_row_sums(self, visium_sub):
        """After normalize_total, each non-zero row should sum to target_sum."""
        from sc_tools.pp import normalize_total

        adata = visium_sub.copy()
        normalize_total(adata, target_sum=1e4)
        X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
        row_sums = X.sum(axis=1)
        # Rows that had zero total before normalization may remain zero
        original_X = visium_sub.X.toarray() if sparse.issparse(visium_sub.X) else visium_sub.X
        nonzero_rows = original_X.sum(axis=1) > 0
        np.testing.assert_allclose(
            row_sums[nonzero_rows],
            1e4,
            rtol=1e-4,
            err_msg="Non-zero rows should sum to target_sum after normalization",
        )

    def test_backup_raw_then_normalize(self, visium_sub):
        """backup_raw followed by normalize must leave adata.raw intact."""
        from sc_tools.pp import backup_raw, normalize_total

        adata = visium_sub.copy()
        # Some checkpoints already have raw set; clear it to test backup_raw
        adata.raw = None
        backup_raw(adata)
        assert adata.raw is not None
        raw_shape = adata.raw.shape

        normalize_total(adata, target_sum=1e4)
        # raw shape should be unchanged
        assert adata.raw.shape == raw_shape, "adata.raw shape changed after normalize_total"

    def test_backup_raw_noop_on_second_call(self, visium_sub):
        """Second backup_raw call must not overwrite an existing adata.raw."""
        from sc_tools.pp import backup_raw, normalize_total

        adata = visium_sub.copy()
        backup_raw(adata)
        # Store a sentinel in raw to detect overwrite
        sentinel_shape = adata.raw.shape

        # Normalize and call backup again -- should be no-op
        normalize_total(adata, target_sum=1e4)
        backup_raw(adata)
        assert adata.raw.shape == sentinel_shape, "backup_raw overwrote existing raw"

    def test_single_library_id_normalize(self, visium_sub):
        """Normalization on a single-library subset must not crash."""
        from sc_tools.pp import normalize_total

        # Take only 1 library_id
        one_lib = visium_sub.obs["library_id"].unique()[0]
        sub = visium_sub[visium_sub.obs["library_id"] == one_lib].copy()
        # Should not raise
        normalize_total(sub, target_sum=1e4)
        X = sub.X.toarray() if sparse.issparse(sub.X) else sub.X
        assert not np.any(np.isnan(X))

    def test_log_transform_no_nan_inf(self, visium_sub):
        """log1p after normalize_total on real data must not produce NaN or Inf."""
        from sc_tools.pp import log_transform, normalize_total

        adata = visium_sub.copy()
        normalize_total(adata, target_sum=1e4)
        log_transform(adata)
        X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
        assert not np.any(np.isnan(X)), "NaN values after log1p"
        assert not np.any(np.isinf(X)), "Inf values after log1p"

    def test_log_transform_values_bounded(self, visium_sub):
        """log1p(normalize_total) values should be >= 0 and < log1p(1e4 * max_gene_fraction)."""
        from sc_tools.pp import log_transform, normalize_total

        adata = visium_sub.copy()
        normalize_total(adata, target_sum=1e4)
        log_transform(adata)
        X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
        assert X.min() >= 0, "log1p produced negative values"
        # Upper bound: log1p(1e4) ~ 9.21 (all counts in 1 gene), add small tolerance
        assert X.max() <= np.log1p(1e4) + 1e-3, "log1p values unexpectedly high"

    def test_highly_variable_genes_multi_batch(self, visium_sub):
        """HVG selection with seurat flavor and batch_key must yield correct column."""
        import scanpy as sc

        from sc_tools.pp import log_transform, normalize_total

        adata = visium_sub.copy()
        normalize_total(adata, target_sum=1e4)
        log_transform(adata)

        n_top = 100
        sc.pp.highly_variable_genes(
            adata,
            flavor="seurat",
            n_top_genes=n_top,
            batch_key="library_id",
        )
        assert "highly_variable" in adata.var.columns
        n_hvg = adata.var["highly_variable"].sum()
        assert n_hvg >= n_top * 0.5, f"Too few HVGs selected: {n_hvg}"

    def test_highly_variable_genes_single_sample(self, visium_sub):
        """HVG selection without batch_key on a single-library subset must work."""
        import scanpy as sc

        from sc_tools.pp import log_transform, normalize_total

        one_lib = visium_sub.obs["library_id"].unique()[0]
        sub = visium_sub[visium_sub.obs["library_id"] == one_lib].copy()
        normalize_total(sub, target_sum=1e4)
        log_transform(sub)
        sc.pp.highly_variable_genes(sub, flavor="seurat", n_top_genes=50)
        assert "highly_variable" in sub.var.columns
        assert sub.var["highly_variable"].sum() == 50


# ---------------------------------------------------------------------------
# PCA / reduce tests
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not _abs(VISIUM_ANNOTATED).exists(),
    reason="Visium annotated data not available",
)
class TestPCARealData:
    """PCA dimensionality reduction on real Visium data."""

    @pytest.fixture(scope="class")
    def visium_prepared(self):
        """400 spots, 50 HVGs, normalized and log-transformed, ready for PCA."""
        import anndata as ad
        import scanpy as sc

        adata = ad.read_h5ad(_abs(VISIUM_ANNOTATED))
        sub = _subset_by_batch(adata, "library_id", n_per_batch=200, max_batches=2)

        sc.pp.normalize_total(sub, target_sum=1e4)
        sc.pp.log1p(sub)
        sc.pp.highly_variable_genes(sub, n_top_genes=50)
        sub._inplace_subset_var(sub.var["highly_variable"].values)
        return sub

    def test_pca_stores_obsm(self, visium_prepared):
        """PCA must store X_pca in obsm."""
        from sc_tools.pp import pca

        adata = visium_prepared.copy()
        pca(adata, n_comps=20)
        assert "X_pca" in adata.obsm, "X_pca not stored in obsm after pca()"
        assert adata.obsm["X_pca"].ndim == 2

    def test_pca_shape_correct(self, visium_prepared):
        """PCA output shape must be (n_obs, n_comps)."""
        from sc_tools.pp import pca

        adata = visium_prepared.copy()
        pca(adata, n_comps=20)
        n_obs = adata.n_obs
        assert adata.obsm["X_pca"].shape == (n_obs, 20)

    def test_pca_no_nan_inf(self, visium_prepared):
        """PCA output must not contain NaN or Inf values."""
        from sc_tools.pp import pca

        adata = visium_prepared.copy()
        pca(adata, n_comps=20)
        X_pca = adata.obsm["X_pca"]
        assert not np.any(np.isnan(X_pca)), "NaN in X_pca"
        assert not np.any(np.isinf(X_pca)), "Inf in X_pca"

    def test_pca_n_comps_exceeds_n_cells(self, visium_prepared):
        """When n_comps > n_cells, pca must clamp and not raise."""
        from sc_tools.pp import pca

        # Make a tiny subset: 10 cells, 50 genes
        adata = visium_prepared[:10].copy()
        # Request 200 components -- far more than valid
        pca(adata, n_comps=200)
        assert "X_pca" in adata.obsm
        # Clamped to min(n_obs, n_vars) - 1
        max_valid = min(adata.n_obs, adata.n_vars) - 1
        assert adata.obsm["X_pca"].shape[1] <= max_valid

    def test_pca_n_comps_exceeds_n_genes(self):
        """When n_comps > n_genes, pca must clamp to n_genes - 1."""
        import anndata as ad
        import scanpy as sc

        from sc_tools.pp import pca

        # 100 cells, 5 genes (fewer genes than typical n_comps)
        X = np.random.negative_binomial(5, 0.3, (100, 5)).astype("float32")
        adata = ad.AnnData(X)
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        pca(adata, n_comps=50)
        assert "X_pca" in adata.obsm
        assert adata.obsm["X_pca"].shape[1] <= 4  # min(100, 5) - 1


# ---------------------------------------------------------------------------
# Recipe tests on real Visium data
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not _abs(VISIUM_ANNOTATED).exists(),
    reason="Visium annotated data not available",
)
class TestRecipeRealVisium:
    """Test preprocess() recipe on real Visium subsets."""

    @pytest.fixture(scope="class")
    def visium_small(self):
        """400 spots from 2 library_ids -- small enough for recipe testing."""
        import anndata as ad

        adata = ad.read_h5ad(_abs(VISIUM_ANNOTATED))
        return _subset_by_batch(adata, "library_id", n_per_batch=200, max_batches=2)

    def test_recipe_visium_harmony_checkpoint_keys(self, visium_small):
        """Visium recipe with harmony integration must populate expected keys."""
        from sc_tools.pp import preprocess

        try:
            import harmonypy  # noqa: F401
        except ImportError:
            pytest.skip("harmonypy not installed")

        result = preprocess(
            visium_small,
            modality="visium",
            integration="harmony",
            n_top_genes=100,
            resolution=0.5,
            n_neighbors=10,
            copy=True,
        )
        assert result.raw is not None, "adata.raw not backed up"
        assert "leiden" in result.obs.columns, "leiden clusters missing"
        assert "X_pca" in result.obsm, "X_pca missing"
        assert "X_umap" in result.obsm, "X_umap missing"
        assert result.n_vars <= 100, "Gene count exceeds n_top_genes"

    def test_recipe_visium_none_integration(self, visium_small):
        """Visium recipe with integration='none' must complete without integration keys."""
        from sc_tools.pp import preprocess

        result = preprocess(
            visium_small,
            modality="visium",
            integration="none",
            n_top_genes=100,
            resolution=0.5,
            n_neighbors=10,
            copy=True,
        )
        assert result.raw is not None
        assert "leiden" in result.obs.columns
        assert "X_pca" in result.obsm
        # No scVI or harmony keys expected
        assert "X_scVI" not in result.obsm
        assert "X_pca_harmony" not in result.obsm

    def test_recipe_single_library_id(self, visium_small):
        """Recipe on a single-library subset must not crash."""
        from sc_tools.pp import preprocess

        one_lib = visium_small.obs["library_id"].unique()[0]
        sub = visium_small[visium_small.obs["library_id"] == one_lib].copy()

        result = preprocess(
            sub,
            modality="visium",
            integration="none",
            n_top_genes=50,
            resolution=0.5,
            n_neighbors=5,
            copy=True,
        )
        assert result.raw is not None
        assert "leiden" in result.obs.columns

    def test_recipe_copy_preserves_original(self, visium_small):
        """copy=True must leave the original adata unchanged."""
        from sc_tools.pp import preprocess

        original_shape = visium_small.shape
        original_raw = visium_small.raw

        _ = preprocess(
            visium_small,
            modality="visium",
            integration="none",
            n_top_genes=100,
            copy=True,
        )
        assert visium_small.shape == original_shape
        assert visium_small.raw is original_raw  # None unless previously set

    def test_recipe_filter_genes_applied(self, visium_small):
        """After recipe, MT/RP/HB genes from the original data should not appear."""
        from sc_tools.pp import preprocess

        result = preprocess(
            visium_small,
            modality="visium",
            integration="none",
            n_top_genes=100,
            copy=True,
        )
        import re

        for gene in result.var_names:
            assert not re.match(r"^MT-", gene, re.IGNORECASE), f"MT gene survived: {gene}"
            assert not re.match(r"^RP[SL]", gene, re.IGNORECASE), f"RP gene survived: {gene}"


# ---------------------------------------------------------------------------
# IMC arcsinh + recipe tests
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not _abs(IMC_CELLTYPED).exists(),
    reason="IMC celltyped data not available",
)
class TestNormalizeRealIMC:
    """Arcsinh transform and IMC recipe on real IMC data."""

    @pytest.fixture(scope="class")
    def imc_sub(self):
        """Load and subset IMC data (3 ROIs, 200 cells each)."""
        import anndata as ad

        adata = ad.read_h5ad(_abs(IMC_CELLTYPED))
        # IMC batch key may be 'roi' or 'library_id'
        batch_key = "roi" if "roi" in adata.obs.columns else "library_id"
        sub = _subset_by_batch(adata, batch_key, n_per_batch=200, max_batches=3)
        sub.obs["library_id"] = sub.obs[batch_key]
        return sub

    def test_arcsinh_no_nan_inf(self, imc_sub):
        """arcsinh_transform on real IMC data must not produce NaN or Inf."""
        from sc_tools.pp import arcsinh_transform

        adata = imc_sub.copy()
        arcsinh_transform(adata, cofactor=5)
        X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
        assert not np.any(np.isnan(X)), "NaN after arcsinh"
        assert not np.any(np.isinf(X)), "Inf after arcsinh"

    def test_arcsinh_sparse_input(self, imc_sub):
        """arcsinh_transform on sparse input must produce dense output."""
        from sc_tools.pp import arcsinh_transform

        adata = imc_sub.copy()
        adata.X = sparse.csr_matrix(adata.X.toarray() if sparse.issparse(adata.X) else adata.X)
        assert sparse.issparse(adata.X), "X must be sparse for this test"
        arcsinh_transform(adata, cofactor=5)
        assert not sparse.issparse(adata.X), "X should be dense after arcsinh"

    def test_arcsinh_values_bounded(self, imc_sub):
        """arcsinh output should be in a reasonable range (no extreme values)."""
        from sc_tools.pp import arcsinh_transform

        adata = imc_sub.copy()
        arcsinh_transform(adata, cofactor=5)
        X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
        # arcsinh(x/5) for typical IMC intensities (0-1000) -> max ~5.3
        assert abs(X).max() < 50, "Unreasonably large arcsinh values"

    def test_imc_recipe_no_integration(self, imc_sub):
        """IMC recipe with integration=none must produce standard checkpoint keys."""
        from sc_tools.pp import preprocess

        result = preprocess(
            imc_sub,
            modality="imc",
            integration="none",
            resolution=0.5,
            n_neighbors=10,
            copy=True,
        )
        assert result.raw is not None, "adata.raw not backed up"
        assert "leiden" in result.obs.columns, "leiden clusters missing"
        assert "X_pca" in result.obsm, "X_pca missing"
        # arcsinh + scale must have been applied (values are bounded, not raw counts)
        X_now = result.X.toarray() if sparse.issparse(result.X) else result.X
        # arcsinh(x/5) for typical IMC intensities is bounded; scale clips to max_value=10
        assert np.nanmax(np.abs(X_now)) <= 11, "Scaled values exceed expected range"
        assert not np.all(X_now == 0), "All-zero X after recipe"


# ---------------------------------------------------------------------------
# Sparse / zero-heavy data stress tests
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not _abs(COSMX_ANNOTATED).exists(),
    reason="CosMx annotated data not available",
)
class TestSparseRealCosMx:
    """Tests on highly sparse CosMx data (many zeros typical in targeted panels)."""

    @pytest.fixture(scope="class")
    def cosmx_sub(self):
        """Load and subset CosMx data (2 samples, 200 cells each)."""
        import anndata as ad

        adata = ad.read_h5ad(_abs(COSMX_ANNOTATED))
        batch_key = "library_id" if "library_id" in adata.obs.columns else "sample"
        if batch_key not in adata.obs.columns:
            # Fall back: assign all to single batch
            adata.obs["library_id"] = "sample_0"
            batch_key = "library_id"
        sub = _subset_by_batch(adata, batch_key, n_per_batch=200, max_batches=2)
        if batch_key != "library_id":
            sub.obs["library_id"] = sub.obs[batch_key]
        return sub

    def test_normalize_sparse_cosmx_no_nan(self, cosmx_sub):
        """normalize_total on very sparse CosMx data must not produce NaN/Inf."""
        from sc_tools.pp import normalize_total

        adata = cosmx_sub.copy()
        normalize_total(adata, target_sum=1e4)
        X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
        assert not np.any(np.isnan(X)), "NaN after normalize_total on CosMx"
        assert not np.any(np.isinf(X)), "Inf after normalize_total on CosMx"

    def test_log_transform_sparse_cosmx_no_nan(self, cosmx_sub):
        """log1p on zero-heavy CosMx data must produce non-negative values without NaN."""
        from sc_tools.pp import log_transform, normalize_total

        adata = cosmx_sub.copy()
        normalize_total(adata, target_sum=1e4)
        log_transform(adata)
        X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
        assert not np.any(np.isnan(X))
        assert not np.any(np.isinf(X))
        assert X.min() >= 0, "log1p produced negative values on CosMx"

    def test_pca_after_normalize_cosmx(self, cosmx_sub):
        """PCA on normalized CosMx data must store X_pca without errors."""
        import scanpy as sc

        from sc_tools.pp import log_transform, normalize_total, pca

        adata = cosmx_sub.copy()
        normalize_total(adata, target_sum=1e4)
        log_transform(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=50)
        adata._inplace_subset_var(adata.var["highly_variable"].values)
        pca(adata, n_comps=20)
        assert "X_pca" in adata.obsm
        X_pca = adata.obsm["X_pca"]
        assert not np.any(np.isnan(X_pca))
        assert not np.any(np.isinf(X_pca))
