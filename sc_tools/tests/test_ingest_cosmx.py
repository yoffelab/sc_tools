"""Unit tests for load_cosmx_sample() in sc_tools.ingest.loaders."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

# ---------------------------------------------------------------------------
# Fixtures: synthetic CosMx-like flat file export (10 cells x 5 genes)
# ---------------------------------------------------------------------------

GENES = ["EPCAM", "CD3E", "CD68", "MKI67", "GAPDH"]
N_CELLS = 10
N_GENES = len(GENES)

RNG = np.random.default_rng(42)


def _make_expression_csv(cosmx_dir: Path) -> Path:
    """Write a minimal AtoMx-style expression CSV (cell x gene)."""
    counts = RNG.integers(0, 50, size=(N_CELLS, N_GENES))
    df = pd.DataFrame(counts, columns=GENES)
    df.insert(0, "cell_id", [f"c{i:04d}" for i in range(N_CELLS)])
    df.insert(1, "fov", [str(i % 3 + 1) for i in range(N_CELLS)])
    path = cosmx_dir / "exprMat_file.csv"
    df.to_csv(path, index=False)
    return path


def _make_metadata_csv(cosmx_dir: Path) -> Path:
    """Write a minimal AtoMx-style cell metadata CSV with coordinates."""
    xs = RNG.uniform(100, 5000, size=N_CELLS).astype(np.float64)
    ys = RNG.uniform(100, 5000, size=N_CELLS).astype(np.float64)
    df = pd.DataFrame(
        {
            "cell_id": [f"c{i:04d}" for i in range(N_CELLS)],
            "fov": [str(i % 3 + 1) for i in range(N_CELLS)],
            "CenterX_local_px": xs,
            "CenterY_local_px": ys,
        }
    )
    path = cosmx_dir / "metadata_file.csv"
    df.to_csv(path, index=False)
    return path


@pytest.fixture()
def cosmx_dir(tmp_path: Path) -> Path:
    """Create a minimal synthetic CosMx AtoMx export directory."""
    _make_expression_csv(tmp_path)
    _make_metadata_csv(tmp_path)
    return tmp_path


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestLoadCosmxSample:
    def test_returns_anndata(self, cosmx_dir):
        from sc_tools.ingest.loaders import load_cosmx_sample

        adata = load_cosmx_sample("sample_A", cosmx_dir, panel_tier="1k")
        assert isinstance(adata, ad.AnnData)

    def test_cell_and_gene_counts(self, cosmx_dir):
        from sc_tools.ingest.loaders import load_cosmx_sample

        adata = load_cosmx_sample("sample_A", cosmx_dir, panel_tier="1k")
        assert adata.n_obs == N_CELLS
        assert adata.n_vars == N_GENES

    def test_obs_sample_column(self, cosmx_dir):
        from sc_tools.ingest.loaders import load_cosmx_sample

        adata = load_cosmx_sample("sample_A", cosmx_dir, panel_tier="1k")
        assert "sample" in adata.obs.columns
        assert (adata.obs["sample"] == "sample_A").all()

    def test_obs_library_id_column(self, cosmx_dir):
        from sc_tools.ingest.loaders import load_cosmx_sample

        adata = load_cosmx_sample("sample_A", cosmx_dir, panel_tier="1k")
        assert "library_id" in adata.obs.columns
        assert (adata.obs["library_id"] == "sample_A").all()

    def test_obs_raw_data_dir_column(self, cosmx_dir):
        from sc_tools.ingest.loaders import load_cosmx_sample

        adata = load_cosmx_sample("sample_A", cosmx_dir, panel_tier="1k")
        assert "raw_data_dir" in adata.obs.columns
        assert (adata.obs["raw_data_dir"] == str(cosmx_dir)).all()

    def test_obsm_spatial_shape(self, cosmx_dir):
        from sc_tools.ingest.loaders import load_cosmx_sample

        adata = load_cosmx_sample("sample_A", cosmx_dir, panel_tier="1k")
        assert "spatial" in adata.obsm
        assert adata.obsm["spatial"].shape == (N_CELLS, 2)

    def test_obsm_spatial_dtype_float32(self, cosmx_dir):
        from sc_tools.ingest.loaders import load_cosmx_sample

        adata = load_cosmx_sample("sample_A", cosmx_dir, panel_tier="1k")
        assert adata.obsm["spatial"].dtype == np.float32

    def test_x_raw_counts_sparse_csr(self, cosmx_dir):
        from sc_tools.ingest.loaders import load_cosmx_sample

        adata = load_cosmx_sample("sample_A", cosmx_dir, panel_tier="1k")
        assert sp.issparse(adata.X)
        assert sp.isspmatrix_csr(adata.X)

    def test_x_values_are_non_negative_integers(self, cosmx_dir):
        from sc_tools.ingest.loaders import load_cosmx_sample

        adata = load_cosmx_sample("sample_A", cosmx_dir, panel_tier="1k")
        x_dense = adata.X.toarray()
        assert (x_dense >= 0).all()
        assert np.issubdtype(adata.X.dtype, np.integer), "X must have an integer dtype"

    def test_gene_names_preserved(self, cosmx_dir):
        from sc_tools.ingest.loaders import load_cosmx_sample

        adata = load_cosmx_sample("sample_A", cosmx_dir, panel_tier="1k")
        assert list(adata.var_names) == GENES

    def test_panel_tier_stored_in_uns(self, cosmx_dir):
        from sc_tools.ingest.loaders import load_cosmx_sample

        adata = load_cosmx_sample("sample_A", cosmx_dir, panel_tier="6k")
        assert adata.uns.get("panel_tier") == "6k"

    def test_missing_cosmx_dir_raises(self, tmp_path):
        from sc_tools.ingest.loaders import load_cosmx_sample

        with pytest.raises(FileNotFoundError, match="cosmx_dir"):
            load_cosmx_sample("sample_A", tmp_path / "nonexistent", panel_tier="1k")

    def test_missing_expression_file_raises(self, tmp_path):
        """Directory exists but no expression matrix found."""
        from sc_tools.ingest.loaders import load_cosmx_sample

        # Create a dir with only the metadata file; no expression matrix
        _make_metadata_csv(tmp_path)
        with pytest.raises(FileNotFoundError, match="expression"):
            load_cosmx_sample("sample_A", tmp_path, panel_tier="1k")

    def test_invalid_panel_tier_raises(self, cosmx_dir):
        from sc_tools.ingest.loaders import load_cosmx_sample

        with pytest.raises(ValueError, match="panel_tier"):
            load_cosmx_sample("sample_A", cosmx_dir, panel_tier="bad_tier")

    def test_exported_from_ingest_package(self):
        from sc_tools.ingest import load_cosmx_sample  # noqa: F401 — import must succeed

    def test_parquet_expression_matrix(self, tmp_path):
        """Loader also accepts Parquet expression files."""
        from sc_tools.ingest.loaders import load_cosmx_sample

        # Write expression as parquet
        counts = RNG.integers(0, 50, size=(N_CELLS, N_GENES))
        df = pd.DataFrame(counts, columns=GENES)
        df.insert(0, "cell_id", [f"c{i:04d}" for i in range(N_CELLS)])
        df.insert(1, "fov", [str(i % 3 + 1) for i in range(N_CELLS)])
        df.to_parquet(tmp_path / "exprMat_file.parquet", index=False)
        _make_metadata_csv(tmp_path)

        adata = load_cosmx_sample("sample_B", tmp_path, panel_tier="full_library")
        assert adata.n_obs == N_CELLS
        assert adata.n_vars == N_GENES

    def test_missing_metadata_spatial_zeros(self, tmp_path):
        """When no metadata file is present, obsm['spatial'] is zero-filled."""
        from sc_tools.ingest.loaders import load_cosmx_sample

        # Only the expression file — no metadata_file.csv
        _make_expression_csv(tmp_path)

        adata = load_cosmx_sample("sample_C", tmp_path, panel_tier="1k")
        assert "spatial" in adata.obsm
        assert adata.obsm["spatial"].shape == (N_CELLS, 2)
        assert np.all(adata.obsm["spatial"] == 0)
