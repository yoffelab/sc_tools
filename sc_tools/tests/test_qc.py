"""
Unit tests for sc_tools.qc (Phase 1).

- metrics: calculate_qc_metrics, filter_cells, filter_genes, highly_variable_genes
- plots: qc_2x2_grid (qc_spatial_multipage requires spatial setup; smoke test only if available)
"""

import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from sc_tools.qc import (
    calculate_qc_metrics,
    filter_cells,
    filter_genes,
    highly_variable_genes,
    plot_highly_variable_genes,
    plot_spatially_variable_genes,
    qc_2x2_grid,
    qc_2x4_pre_post,
    qc_scatter_counts_genes,
    qc_violin_metrics,
)


def _minimal_adata(n_obs=80, n_vars=100, mt_genes=5):
    """Synthetic adata with MT- genes for QC metrics."""
    np.random.seed(42)
    X = np.random.negative_binomial(5, 0.3, (n_obs, n_vars)).astype(np.float32)
    var_names = [f"MT-{i}" for i in range(mt_genes)] + [f"g{i}" for i in range(mt_genes, n_vars)]
    adata = sc.AnnData(
        X,
        obs=pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)]),
        var=pd.DataFrame(index=var_names),
    )
    return adata


def test_calculate_qc_metrics():
    adata = _minimal_adata()
    calculate_qc_metrics(adata, inplace=True, percent_top=(10, 20, 50))
    assert "total_counts" in adata.obs.columns
    assert "n_genes_by_counts" in adata.obs.columns
    assert "pct_counts_mt" in adata.obs.columns
    assert "mt" in adata.var.columns
    assert adata.var["mt"].sum() == 5


def test_filter_cells_genes():
    adata = _minimal_adata()
    calculate_qc_metrics(adata, inplace=True, percent_top=(10, 20))
    n_obs, n_vars = adata.n_obs, adata.n_vars
    filter_cells(adata, min_genes=3, inplace=True)
    filter_genes(adata, min_cells=3, inplace=True)
    assert adata.n_obs <= n_obs
    assert adata.n_vars <= n_vars


def test_highly_variable_genes():
    adata = _minimal_adata()
    calculate_qc_metrics(adata, inplace=True, percent_top=(10, 20))
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    highly_variable_genes(adata, flavor="seurat", n_top_genes=20, inplace=True)
    assert "highly_variable" in adata.var.columns
    assert adata.var["highly_variable"].sum() >= 20  # may exceed slightly due to ties


def test_qc_2x2_grid():
    adata = _minimal_adata()
    calculate_qc_metrics(adata, inplace=True, percent_top=(10, 20))
    fig = qc_2x2_grid(adata)
    assert fig is not None
    assert len(fig.axes) == 4
    fig.clear()
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_qc_2x2_grid_save(tmp_path):
    adata = _minimal_adata()
    calculate_qc_metrics(adata, inplace=True, percent_top=(10, 20))
    qc_2x2_grid(adata, output_dir=tmp_path, basename="qc_test")
    assert (tmp_path / "qc_test.pdf").exists()
    assert (tmp_path / "qc_test.png").exists()


def test_qc_2x4_pre_post(tmp_path):
    adata_pre = _minimal_adata(n_obs=100, n_vars=80)
    calculate_qc_metrics(adata_pre, inplace=True, percent_top=(10, 20))
    adata_post = _minimal_adata(n_obs=60, n_vars=70)  # fewer obs/vars as if filtered
    calculate_qc_metrics(adata_post, inplace=True, percent_top=(10, 20))
    fig = qc_2x4_pre_post(adata_pre, adata_post, output_dir=tmp_path, basename="qc_pre_post")
    assert fig is not None
    assert len(fig.axes) == 8
    assert (tmp_path / "qc_pre_post.pdf").exists()
    assert (tmp_path / "qc_pre_post.png").exists()
    fig.clear()
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_qc_violin_metrics(tmp_path):
    adata = _minimal_adata()
    calculate_qc_metrics(adata, inplace=True, percent_top=(10, 20))
    fig = qc_violin_metrics(adata, output_dir=tmp_path, basename="violin")
    assert fig is not None
    qc_violin_metrics(adata, output_dir=tmp_path, basename="violin2")
    assert (tmp_path / "violin2.pdf").exists()
    assert (tmp_path / "violin2.png").exists()
    fig.clear()
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_qc_scatter_counts_genes(tmp_path):
    adata = _minimal_adata()
    calculate_qc_metrics(adata, inplace=True, percent_top=(10, 20))
    fig = qc_scatter_counts_genes(adata, output_dir=tmp_path, basename="scatter")
    assert fig is not None
    assert (tmp_path / "scatter.pdf").exists()
    fig.clear()
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_plot_highly_variable_genes(tmp_path):
    adata = _minimal_adata()
    calculate_qc_metrics(adata, inplace=True, percent_top=(10, 20))
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    highly_variable_genes(adata, flavor="seurat", n_top_genes=20, inplace=True)
    fig = plot_highly_variable_genes(adata, output_dir=tmp_path, basename="hvg")
    assert fig is not None
    assert (tmp_path / "hvg.pdf").exists()
    fig.clear()
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_plot_spatially_variable_genes_no_data(tmp_path):
    """plot_spatially_variable_genes with no spatial_i returns a placeholder figure."""
    adata = _minimal_adata()
    fig = plot_spatially_variable_genes(adata, output_dir=tmp_path, basename="svg")
    assert fig is not None
    assert (tmp_path / "svg.pdf").exists()
    fig.clear()
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_plot_spatially_variable_genes_with_var(tmp_path):
    """plot_spatially_variable_genes with spatial_i in var (no squidpy)."""
    adata = _minimal_adata(n_obs=50, n_vars=30)
    adata.var["spatial_i"] = np.random.RandomState(42).rand(adata.n_vars).astype(float)
    adata.var["spatial_pval"] = np.random.RandomState(43).rand(adata.n_vars).astype(float)
    adata.var["spatially_variable"] = adata.var["spatial_i"] > 0.5
    fig = plot_spatially_variable_genes(adata, output_dir=tmp_path, basename="svg_var")
    assert fig is not None
    assert (tmp_path / "svg_var.pdf").exists()
    fig.clear()
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_spatially_variable_genes_per_library_no_library_id():
    """When library_id_col is missing from obs, returns None without error."""
    from sc_tools.qc import spatially_variable_genes_per_library

    adata = _minimal_adata(n_obs=20, n_vars=15)
    assert "library_id" not in adata.obs.columns
    out = spatially_variable_genes_per_library(adata, library_id_col="library_id")
    assert out is None


@pytest.mark.skip(reason="squidpy optional; run when squidpy installed")
def test_spatially_variable_genes():
    from sc_tools.qc import spatially_variable_genes

    adata = _minimal_adata(n_obs=60, n_vars=30)
    adata.obsm["spatial"] = np.column_stack(
        [np.repeat(np.arange(6), 10), np.tile(np.arange(10), 6)]
    )
    spatially_variable_genes(
        adata, genes=adata.var_names[:5].tolist(), n_top_genes=2, copy_var=True
    )
    assert "spatial_i" in adata.var.columns or "spatially_variable" in adata.var.columns
