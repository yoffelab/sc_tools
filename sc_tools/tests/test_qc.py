"""
Unit tests for sc_tools.qc (Phase 1).

- metrics: calculate_qc_metrics, filter_cells, filter_genes, highly_variable_genes
- plots: qc_2x2_grid, cross-sample comparison plots
- sample_qc: filter_spots, compute_sample_metrics, classify_samples, apply_qc_filter
- report: generate_qc_report
"""

import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from sc_tools.qc import (
    apply_qc_filter,
    calculate_qc_metrics,
    classify_samples,
    compute_sample_metrics,
    filter_cells,
    filter_genes,
    filter_spots,
    highly_variable_genes,
    plot_highly_variable_genes,
    plot_spatially_variable_genes,
    qc_2x2_grid,
    qc_2x4_pre_post,
    qc_sample_comparison_bar,
    qc_sample_scatter_matrix,
    qc_sample_violin_grouped,
    qc_scatter_counts_genes,
    qc_violin_metrics,
    save_pass_fail_lists,
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


def _multi_sample_adata(n_samples=6, n_obs_per=50, n_vars=100, mt_genes=5):
    """Synthetic multi-sample adata for sample-level QC."""
    np.random.seed(42)
    adatas = []
    for i in range(n_samples):
        X = np.random.negative_binomial(5, 0.3, (n_obs_per, n_vars)).astype(np.float32)
        var_names = [f"MT-{j}" for j in range(mt_genes)] + [
            f"g{j}" for j in range(mt_genes, n_vars)
        ]
        obs = pd.DataFrame(
            {"library_id": f"sample_{i}"},
            index=[f"cell_{i}_{j}" for j in range(n_obs_per)],
        )
        ad = sc.AnnData(X, obs=obs, var=pd.DataFrame(index=var_names))
        adatas.append(ad)
    import anndata

    adata = anndata.concat(adatas, join="outer")
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    return adata


# ---------------------------------------------------------------------------
# Existing metrics tests
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# Existing plot tests
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# filter_spots tests
# ---------------------------------------------------------------------------


def test_filter_spots_visium_hd():
    """Spot filtering with Visium HD defaults (min_counts=10, min_genes=5)."""
    adata = _minimal_adata(n_obs=100, n_vars=50)
    # Artificially set some spots to very low counts
    adata.X[:10, :] = 0  # 10 spots with zero counts
    adata.X[:10, :3] = 1  # but 3 genes with 1 count each (total_counts=3, n_genes=3)
    calculate_qc_metrics(adata, inplace=True, percent_top=None)
    n_before = adata.n_obs
    filter_spots(adata, modality="visium_hd", inplace=True)
    assert adata.n_obs < n_before
    # All remaining spots should pass thresholds
    assert (adata.obs["total_counts"] >= 10).all() or (adata.obs["n_genes_by_counts"] >= 5).all()


def test_filter_spots_preserves_good():
    """Good spots are not removed by filter_spots."""
    adata = _minimal_adata(n_obs=50, n_vars=80)
    calculate_qc_metrics(adata, inplace=True, percent_top=None)
    # With negative binomial(5, 0.3), counts should be high enough
    n_before = adata.n_obs
    filter_spots(adata, modality="visium", inplace=True)
    # Most/all spots should survive
    assert adata.n_obs == n_before  # negative_binomial(5,0.3) gives high counts


def test_filter_spots_copy():
    """filter_spots with inplace=False returns a copy."""
    adata = _minimal_adata(n_obs=50, n_vars=30)
    adata.X[:5, :] = 0
    calculate_qc_metrics(adata, inplace=True, percent_top=None)
    result = filter_spots(adata, modality="visium", inplace=False)
    assert result is not None
    assert result.n_obs <= adata.n_obs


# ---------------------------------------------------------------------------
# compute_sample_metrics tests
# ---------------------------------------------------------------------------


def test_compute_sample_metrics():
    """Compute per-sample metrics from multi-sample AnnData."""
    adata = _multi_sample_adata(n_samples=4, n_obs_per=30)
    calculate_qc_metrics(adata, inplace=True, percent_top=None)
    metrics = compute_sample_metrics(adata, sample_col="library_id")
    assert len(metrics) == 4
    assert "n_spots" in metrics.columns
    assert "total_counts_median" in metrics.columns
    assert "n_genes_median" in metrics.columns
    assert "n_genes_detected" in metrics.columns
    assert "pct_mt_median" in metrics.columns
    assert "pct_mt_gt5" in metrics.columns
    # All samples should have 30 spots
    assert (metrics["n_spots"] == 30).all()


def test_compute_sample_metrics_missing_col():
    """Raises ValueError when sample_col is not in obs."""
    adata = _minimal_adata()
    with pytest.raises(ValueError, match="sample_col"):
        compute_sample_metrics(adata, sample_col="nonexistent")


# ---------------------------------------------------------------------------
# classify_samples tests
# ---------------------------------------------------------------------------


def test_classify_samples_absolute():
    """Samples below absolute thresholds are flagged."""
    metrics = pd.DataFrame(
        {
            "n_spots": [100, 10, 200, 150],
            "n_genes_median": [500, 5, 300, 400],
            "total_counts_median": [2000, 10, 1500, 1800],
            "pct_mt_median": [3.0, 60.0, 5.0, 4.0],
        },
        index=["s1", "s2", "s3", "s4"],
    )
    classified = classify_samples(metrics, modality="visium")
    assert not classified.loc["s2", "qc_pass"]
    assert classified.loc["s2", "qc_flag_absolute"]
    assert classified.loc["s1", "qc_pass"]
    assert classified.loc["s3", "qc_pass"]


def test_classify_samples_outlier():
    """MAD-based outlier detection flags extreme samples."""
    np.random.seed(42)
    n = 25
    metrics = pd.DataFrame(
        {
            "n_spots": np.concatenate([[500] * (n - 1), [10]]),
            "n_genes_median": np.concatenate([[300] * (n - 1), [5]]),
            "total_counts_median": np.concatenate([[2000] * (n - 1), [50]]),
            "pct_mt_median": np.concatenate([[3.0] * (n - 1), [3.0]]),
        },
        index=[f"s{i}" for i in range(n)],
    )
    classified = classify_samples(metrics, modality="visium", mad_multiplier=3.0)
    # The last sample should be flagged as outlier (extreme low values)
    # It may also be flagged by absolute thresholds
    assert not classified.loc[f"s{n - 1}", "qc_pass"]


def test_classify_samples_small_cohort():
    """With <10 samples, outlier detection is very conservative (mad_multiplier=5.0)."""
    metrics = pd.DataFrame(
        {
            "n_spots": [200, 150, 180, 50],
            "n_genes_median": [300, 250, 280, 100],
            "total_counts_median": [2000, 1800, 1900, 500],
            "pct_mt_median": [3.0, 4.0, 3.5, 10.0],
        },
        index=["s1", "s2", "s3", "s4"],
    )
    classified = classify_samples(metrics, modality="visium")
    # With only 4 samples, outlier detection should not fire (MAD mult=5.0)
    # Absolute thresholds: s4 has n_genes_median=100 > 50 min, so passes absolute
    assert classified.loc["s1", "qc_pass"]
    assert classified.loc["s2", "qc_pass"]
    # s4 should pass absolute thresholds (100 > 50 min genes)
    # and outlier detection should be skipped (n<5 by default)
    assert classified.loc["s4", "qc_pass"]


def test_classify_samples_all_pass():
    """Normal cohort where all samples pass."""
    metrics = pd.DataFrame(
        {
            "n_spots": [500, 450, 480, 520, 490, 510],
            "n_genes_median": [300, 280, 310, 290, 305, 295],
            "total_counts_median": [2000, 1900, 2100, 2050, 1950, 2000],
            "pct_mt_median": [3.0, 3.5, 2.8, 3.2, 3.1, 3.3],
        },
        index=[f"s{i}" for i in range(6)],
    )
    classified = classify_samples(metrics, modality="visium")
    assert classified["qc_pass"].all()
    assert (classified["qc_fail_reasons"] == "").all()


# ---------------------------------------------------------------------------
# save_pass_fail_lists tests
# ---------------------------------------------------------------------------


def test_save_pass_fail_lists(tmp_path):
    """CSV files are created with correct content."""
    metrics = pd.DataFrame(
        {
            "n_spots": [200, 10],
            "n_genes_median": [300, 5],
            "total_counts_median": [2000, 10],
            "pct_mt_median": [3.0, 60.0],
        },
        index=["good", "bad"],
    )
    classified = classify_samples(metrics, modality="visium")
    pass_path, fail_path = save_pass_fail_lists(classified, tmp_path)
    assert pass_path.exists()
    assert fail_path.exists()
    pass_df = pd.read_csv(pass_path, index_col=0)
    fail_df = pd.read_csv(fail_path, index_col=0)
    assert "good" in pass_df.index
    assert "bad" in fail_df.index


# ---------------------------------------------------------------------------
# apply_qc_filter tests
# ---------------------------------------------------------------------------


def test_apply_qc_filter_backup(tmp_path):
    """Backup file is created and filtered file has fewer obs."""
    adata = _multi_sample_adata(n_samples=3, n_obs_per=40)
    calculate_qc_metrics(adata, inplace=True, percent_top=None)
    metrics = compute_sample_metrics(adata, sample_col="library_id")

    # Manually set one sample to fail
    classified = classify_samples(metrics, modality="visium")
    classified.loc["sample_0", "qc_pass"] = False
    classified.loc["sample_0", "qc_flag_absolute"] = True
    classified.loc["sample_0", "qc_fail_reasons"] = "manually failed for test"

    backup = tmp_path / "backup.h5ad"
    output = tmp_path / "filtered.h5ad"

    result = apply_qc_filter(
        adata,
        classified,
        sample_col="library_id",
        modality="visium",
        output_path=output,
        backup_path=backup,
    )
    assert backup.exists()
    assert output.exists()
    # sample_0 removed: should have 2 * 40 = 80 spots (minus any spot filtering)
    assert result.n_obs <= 80
    assert "sample_0" not in result.obs["library_id"].values


# ---------------------------------------------------------------------------
# Cross-sample comparison plot tests
# ---------------------------------------------------------------------------


def test_qc_sample_comparison_bar(tmp_path):
    """Bar chart generates without error."""
    import matplotlib.pyplot as plt

    metrics = pd.DataFrame(
        {
            "n_spots": [200, 150, 300],
            "n_genes_median": [300, 250, 350],
            "total_counts_median": [2000, 1500, 2500],
        },
        index=["s1", "s2", "s3"],
    )
    fig = qc_sample_comparison_bar(metrics, output_dir=tmp_path, basename="bar_test")
    assert fig is not None
    assert (tmp_path / "bar_test.png").exists()
    plt.close(fig)


def test_qc_sample_violin_grouped(tmp_path):
    """Violin grouped by sample generates without error."""
    import matplotlib.pyplot as plt

    adata = _multi_sample_adata(n_samples=3, n_obs_per=30)
    calculate_qc_metrics(adata, inplace=True, percent_top=None)
    fig = qc_sample_violin_grouped(
        adata,
        sample_col="library_id",
        output_dir=tmp_path,
        basename="violin_grouped",
    )
    assert fig is not None
    assert (tmp_path / "violin_grouped.png").exists()
    plt.close(fig)


def test_qc_sample_scatter_matrix(tmp_path):
    """Scatter matrix generates without error."""
    import matplotlib.pyplot as plt

    metrics = pd.DataFrame(
        {
            "n_spots": [200, 150, 300, 250],
            "n_genes_median": [300, 250, 350, 280],
            "total_counts_median": [2000, 1500, 2500, 2200],
            "pct_mt_median": [3.0, 4.0, 2.5, 3.5],
        },
        index=["s1", "s2", "s3", "s4"],
    )
    fig = qc_sample_scatter_matrix(metrics, output_dir=tmp_path, basename="scatter_matrix")
    assert fig is not None
    assert (tmp_path / "scatter_matrix.png").exists()
    plt.close(fig)


# ---------------------------------------------------------------------------
# generate_qc_report test
# ---------------------------------------------------------------------------


def test_generate_qc_report(tmp_path):
    """HTML report is generated successfully."""
    pytest.importorskip("jinja2")
    from sc_tools.qc import generate_qc_report

    adata = _multi_sample_adata(n_samples=3, n_obs_per=30)
    calculate_qc_metrics(adata, inplace=True, percent_top=None)
    metrics = compute_sample_metrics(adata, sample_col="library_id")
    classified = classify_samples(metrics, modality="visium")

    # Create a fake figures dir with a QC subfolder
    figures_dir = tmp_path / "figures"
    qc_dir = figures_dir / "QC"
    qc_dir.mkdir(parents=True)

    output = qc_dir / "qc_report.html"
    result = generate_qc_report(
        adata,
        metrics,
        classified,
        figures_dir,
        output,
        sample_col="library_id",
        modality="visium",
    )
    assert result.exists()
    content = result.read_text()
    assert "QC Report" in content
    assert "sample_0" in content
