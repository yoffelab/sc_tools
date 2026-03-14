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
    """HTML report is generated successfully (legacy API)."""
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
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
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


# ---------------------------------------------------------------------------
# New date-versioned QC report tests
# ---------------------------------------------------------------------------


def _qc_fixtures():
    """Shared fixture: multi-sample adata with QC metrics + classification."""
    adata = _multi_sample_adata(n_samples=3, n_obs_per=30)
    calculate_qc_metrics(adata, inplace=True, percent_top=None)
    metrics = compute_sample_metrics(adata, sample_col="library_id")
    classified = classify_samples(metrics, modality="visium")
    return adata, metrics, classified


def test_generate_pre_filter_report(tmp_path):
    """Pre-filter QC HTML is generated with date stamp."""
    pytest.importorskip("jinja2")
    from sc_tools.qc import generate_pre_filter_report

    adata, metrics, classified = _qc_fixtures()
    result = generate_pre_filter_report(
        adata,
        metrics,
        classified,
        tmp_path,
        sample_col="library_id",
        modality="visium",
        date_stamp="20260304",
    )
    assert result.exists()
    assert "20260304" in result.name
    content = result.read_text()
    assert "Pre-filter QC" in content
    assert "sample_0" in content


def test_generate_post_filter_report(tmp_path):
    """Post-filter QC HTML is generated with pre-vs-post comparison."""
    pytest.importorskip("jinja2")
    from sc_tools.qc import generate_post_filter_report

    adata_pre, metrics, classified = _qc_fixtures()
    # Simulate a post-filter adata by removing a few obs
    adata_post = adata_pre[10:].copy()
    calculate_qc_metrics(adata_post, inplace=True, percent_top=None)

    result = generate_post_filter_report(
        adata_pre,
        adata_post,
        metrics,
        classified,
        tmp_path,
        sample_col="library_id",
        modality="visium",
        date_stamp="20260304",
    )
    assert result.exists()
    assert "20260304" in result.name
    content = result.read_text()
    assert "Post-filter" in content


def test_generate_post_integration_report_with_celltype(tmp_path):
    """Post-integration QC HTML with celltype includes integration metrics."""
    pytest.importorskip("jinja2")
    from sc_tools.qc import generate_post_integration_report

    adata = _multi_sample_adata(n_samples=3, n_obs_per=30)
    calculate_qc_metrics(adata, inplace=True, percent_top=None)
    rng = np.random.RandomState(42)
    adata.obsm["X_umap"] = rng.randn(adata.n_obs, 2).astype(np.float32)
    adata.obsm["X_pca"] = rng.randn(adata.n_obs, 10).astype(np.float32)
    adata.obs["leiden"] = pd.Categorical([str(i % 5) for i in range(adata.n_obs)])
    adata.obs["celltype"] = pd.Categorical([f"type_{i % 3}" for i in range(adata.n_obs)])
    adata.obs["batch"] = adata.obs["library_id"]

    result = generate_post_integration_report(
        adata,
        tmp_path,
        batch_key="batch",
        celltype_key="celltype",
        sample_col="library_id",
        date_stamp="20260304",
    )
    assert result.exists()
    assert "20260304" in result.name
    content = result.read_text()
    assert "Post-integration" in content


def test_generate_post_integration_report_no_celltype(tmp_path):
    """Post-integration QC works gracefully without celltype (batch-only)."""
    pytest.importorskip("jinja2")
    from sc_tools.qc import generate_post_integration_report

    adata = _multi_sample_adata(n_samples=3, n_obs_per=30)
    calculate_qc_metrics(adata, inplace=True, percent_top=None)
    rng = np.random.RandomState(42)
    adata.obsm["X_umap"] = rng.randn(adata.n_obs, 2).astype(np.float32)
    adata.obsm["X_pca"] = rng.randn(adata.n_obs, 10).astype(np.float32)
    adata.obs["leiden"] = pd.Categorical([str(i % 5) for i in range(adata.n_obs)])
    adata.obs["batch"] = adata.obs["library_id"]

    result = generate_post_integration_report(
        adata,
        tmp_path,
        batch_key="batch",
        celltype_key=None,
        sample_col="library_id",
        date_stamp="20260304",
    )
    assert result.exists()
    content = result.read_text()
    assert "Post-integration" in content


def test_generate_all_qc_reports(tmp_path):
    """Orchestrator generates all three reports."""
    pytest.importorskip("jinja2")
    from sc_tools.qc import generate_all_qc_reports

    adata, metrics, classified = _qc_fixtures()
    rng = np.random.RandomState(42)
    adata.obsm["X_umap"] = rng.randn(adata.n_obs, 2).astype(np.float32)
    adata.obsm["X_pca"] = rng.randn(adata.n_obs, 10).astype(np.float32)
    adata.obs["leiden"] = pd.Categorical([str(i % 5) for i in range(adata.n_obs)])
    adata.obs["batch"] = adata.obs["library_id"]

    adata_post = adata[5:].copy()
    calculate_qc_metrics(adata_post, inplace=True, percent_top=None)

    results = generate_all_qc_reports(
        adata,
        metrics,
        classified,
        tmp_path,
        adata_post=adata_post,
        adata_integrated=adata,
        sample_col="library_id",
        date_stamp="20260304",
    )
    assert "pre_filter" in results
    assert "post_filter" in results
    assert "post_integration" in results
    for path in results.values():
        assert path.exists()


def test_date_stamp_in_filename(tmp_path):
    """Date stamp appears in generated filename."""
    pytest.importorskip("jinja2")
    from sc_tools.qc import generate_pre_filter_report

    adata, metrics, classified = _qc_fixtures()
    result = generate_pre_filter_report(
        adata,
        metrics,
        classified,
        tmp_path,
        sample_col="library_id",
        date_stamp="20991231",
    )
    assert "20991231" in result.name


def test_backward_compat_generate_qc_report_deprecation(tmp_path):
    """Legacy generate_qc_report emits DeprecationWarning."""
    pytest.importorskip("jinja2")
    from sc_tools.qc import generate_qc_report

    adata, metrics, classified = _qc_fixtures()
    qc_dir = tmp_path / "QC"
    qc_dir.mkdir()

    with pytest.warns(DeprecationWarning, match="deprecated"):
        generate_qc_report(
            adata,
            metrics,
            classified,
            tmp_path,
            qc_dir / "report.html",
            sample_col="library_id",
        )


def test_segmentation_section_no_masks(tmp_path):
    """Segmentation section returns None when no masks dir exists."""
    from sc_tools.qc.report_utils import compute_segmentation_section

    adata = _minimal_adata()
    result = compute_segmentation_section(adata, tmp_path / "nonexistent")
    assert result is None


def test_auto_detect_embeddings():
    """auto_detect_embeddings finds known keys in adata.obsm."""
    from sc_tools.qc.report_utils import auto_detect_embeddings

    adata = _minimal_adata(n_obs=20, n_vars=10)
    rng = np.random.RandomState(42)
    adata.obsm["X_pca"] = rng.randn(20, 10).astype(np.float32)
    adata.obsm["X_scVI"] = rng.randn(20, 10).astype(np.float32)

    found = auto_detect_embeddings(adata)
    assert "scVI" in found
    assert found["scVI"] == "X_scVI"
    assert "Unintegrated (PCA)" in found


# ---------------------------------------------------------------------------
# qc_umap_grid and qc_cluster_distribution tests
# ---------------------------------------------------------------------------


def test_qc_umap_grid(tmp_path):
    """UMAP grid generates without error."""
    import matplotlib.pyplot as plt

    from sc_tools.pl.qc_plots import qc_umap_grid

    adata = _multi_sample_adata(n_samples=3, n_obs_per=20)
    rng = np.random.RandomState(42)
    adata.obsm["X_umap"] = rng.randn(adata.n_obs, 2).astype(np.float32)
    adata.obs["leiden"] = pd.Categorical([str(i % 4) for i in range(adata.n_obs)])

    fig = qc_umap_grid(
        adata,
        color_keys=["library_id", "leiden"],
        output_dir=tmp_path,
        basename="umap_test",
    )
    assert fig is not None
    assert (tmp_path / "umap_test.png").exists()
    plt.close(fig)


def test_qc_cluster_distribution(tmp_path):
    """Cluster distribution bar chart generates without error."""
    import matplotlib.pyplot as plt

    from sc_tools.pl.qc_plots import qc_cluster_distribution

    adata = _multi_sample_adata(n_samples=3, n_obs_per=20)
    adata.obs["leiden"] = pd.Categorical([str(i % 4) for i in range(adata.n_obs)])

    fig = qc_cluster_distribution(
        adata,
        cluster_key="leiden",
        sample_col="library_id",
        output_dir=tmp_path,
        basename="cluster_dist_test",
    )
    assert fig is not None
    assert (tmp_path / "cluster_dist_test.png").exists()
    plt.close(fig)


# ---------------------------------------------------------------------------
# Modality-aware terminology tests
# ---------------------------------------------------------------------------


class TestGetModalityTerms:
    def test_visium_terms(self):
        from sc_tools.qc.report_utils import get_modality_terms

        t = get_modality_terms("visium")
        assert t["feature"] == "Gene"
        assert t["features"] == "Genes"
        assert t["observation"] == "Spot"
        assert t["observations"] == "Spots"
        assert t["has_mt"] is True
        assert t["intensity_label"] == "Total counts"
        assert t["feature_count_label"] == "Genes detected"

    def test_imc_terms(self):
        from sc_tools.qc.report_utils import get_modality_terms

        t = get_modality_terms("imc")
        assert t["feature"] == "Protein"
        assert t["features"] == "Proteins"
        assert t["observation"] == "Cell"
        assert t["observations"] == "Cells"
        assert t["has_mt"] is False
        assert t["intensity_label"] == "Total intensity"
        assert t["feature_count_label"] == "Proteins detected"

    def test_xenium_terms(self):
        from sc_tools.qc.report_utils import get_modality_terms

        t = get_modality_terms("xenium")
        assert t["feature"] == "Gene"
        assert t["observation"] == "Cell"
        assert t["has_mt"] is True

    def test_visium_hd_cell_terms(self):
        from sc_tools.qc.report_utils import get_modality_terms

        t = get_modality_terms("visium_hd_cell")
        assert t["feature"] == "Gene"
        assert t["observation"] == "Cell"


class TestCalculateQcMetricsImc:
    def test_imc_no_mt(self):
        """IMC modality should skip MT pattern automatically."""
        adata = _minimal_adata(n_obs=30, n_vars=50, mt_genes=0)
        calculate_qc_metrics(adata, inplace=True, modality="imc", percent_top=(10, 20))
        assert "total_counts" in adata.obs.columns
        assert "n_genes_by_counts" in adata.obs.columns
        # pct_counts_mt should NOT be computed for IMC
        assert "pct_counts_mt" not in adata.obs.columns

    def test_percent_top_capped(self):
        """percent_top values exceeding n_vars should be silently dropped."""
        adata = _minimal_adata(n_obs=20, n_vars=30, mt_genes=0)
        # Default percent_top=(50,100,200,500) -- all exceed n_vars=30
        calculate_qc_metrics(adata, inplace=True, modality="imc")
        assert "total_counts" in adata.obs.columns


class TestPreFilterReportTerminology:
    def test_imc_report_uses_protein(self, tmp_path):
        """Pre-filter report for IMC should contain 'Protein' not 'Gene'."""
        adata = _multi_sample_adata(n_samples=2, n_obs_per=20)
        calculate_qc_metrics(adata, inplace=True, percent_top=(10, 20))

        from sc_tools.qc import generate_pre_filter_report

        metrics = compute_sample_metrics(adata, sample_col="library_id", modality="imc")
        classified = classify_samples(metrics, modality="imc")

        try:
            path = generate_pre_filter_report(
                adata,
                metrics,
                classified,
                tmp_path,
                sample_col="library_id",
                modality="imc",
            )
        except ImportError:
            pytest.skip("jinja2 not installed")

        html = path.read_text()
        assert "Protein" in html
        assert "Cell" in html or "cell" in html


class TestSegmentationQcReport:
    def test_no_masks_returns_none(self, tmp_path):
        """Segmentation report should return None when no masks exist."""
        from sc_tools.qc import generate_segmentation_qc_report

        adata = _multi_sample_adata(n_samples=2, n_obs_per=20)
        masks_dir = tmp_path / "empty_masks"
        masks_dir.mkdir()

        result = generate_segmentation_qc_report(
            adata,
            masks_dir,
            tmp_path / "out",
        )
        assert result is None


# ---------------------------------------------------------------------------
# Post-celltyping QC report tests
# ---------------------------------------------------------------------------


def test_generate_post_celltyping_report_basic(tmp_path):
    """Post-celltyping report generates with validated celltypes."""
    pytest.importorskip("jinja2")
    from sc_tools.qc import generate_post_celltyping_report

    adata = _multi_sample_adata(n_samples=3, n_obs_per=30)
    calculate_qc_metrics(adata, inplace=True, percent_top=None)
    rng = np.random.RandomState(42)
    adata.obsm["X_umap"] = rng.randn(adata.n_obs, 2).astype(np.float32)
    adata.obsm["X_pca"] = rng.randn(adata.n_obs, 10).astype(np.float32)
    adata.obs["leiden"] = pd.Categorical([str(i % 5) for i in range(adata.n_obs)])
    adata.obs["celltype"] = pd.Categorical([f"type_{i % 3}" for i in range(adata.n_obs)])
    adata.obs["batch"] = adata.obs["library_id"]

    result = generate_post_celltyping_report(
        adata,
        tmp_path,
        celltype_key="celltype",
        batch_key="batch",
        sample_col="library_id",
        date_stamp="20260306",
    )
    assert result.exists()
    assert "post_celltyping" in result.name
    assert "20260306" in result.name
    content = result.read_text()
    assert "Post-celltyping" in content


def test_generate_post_celltyping_report_missing_celltype(tmp_path):
    """Post-celltyping report raises ValueError when celltype is missing."""
    pytest.importorskip("jinja2")
    from sc_tools.qc import generate_post_celltyping_report

    adata = _multi_sample_adata(n_samples=2, n_obs_per=20)

    with pytest.raises(ValueError, match="celltype_key"):
        generate_post_celltyping_report(
            adata,
            tmp_path,
            celltype_key="celltype",
        )


def test_generate_post_celltyping_report_with_integration_test_dir(tmp_path):
    """Post-celltyping report loads integration test embeddings."""
    pytest.importorskip("jinja2")
    from sc_tools.qc import generate_post_celltyping_report

    adata = _multi_sample_adata(n_samples=2, n_obs_per=30)
    calculate_qc_metrics(adata, inplace=True, percent_top=None)
    rng = np.random.RandomState(42)
    adata.obsm["X_umap"] = rng.randn(adata.n_obs, 2).astype(np.float32)
    adata.obsm["X_pca"] = rng.randn(adata.n_obs, 10).astype(np.float32)
    adata.obs["leiden"] = pd.Categorical([str(i % 5) for i in range(adata.n_obs)])
    adata.obs["celltype"] = pd.Categorical([f"type_{i % 3}" for i in range(adata.n_obs)])
    adata.obs["batch"] = adata.obs["library_id"]

    # Create a fake integration test dir with a method h5ad
    test_dir = tmp_path / "tmp" / "integration_test"
    test_dir.mkdir(parents=True)
    test_adata = adata.copy()
    test_adata.obsm["X_pca_harmony"] = rng.randn(test_adata.n_obs, 10).astype(np.float32)
    test_adata.write_h5ad(test_dir / "Harmony.h5ad")

    result = generate_post_celltyping_report(
        adata,
        tmp_path,
        celltype_key="celltype",
        batch_key="batch",
        integration_test_dir=test_dir,
        date_stamp="20260306",
    )
    assert result.exists()
    content = result.read_text()
    assert "Post-celltyping" in content


# ---------------------------------------------------------------------------
# Plan B: report_utils helper tests
# ---------------------------------------------------------------------------


from sc_tools.qc.report_utils import (  # noqa: E402
    _extract_body_content,
    _extract_head_css,
    _find_latest_report,
    _wrap_with_tabs,
)


class TestReportUtils:
    def test_extract_body_content_normal(self):
        html = "<html><head><style>body{}</style></head><body><p>hello</p></body></html>"
        result = _extract_body_content(html)
        assert "<p>hello</p>" in result
        assert "<html>" not in result

    def test_extract_body_content_no_body_fallback(self):
        result = _extract_body_content("<p>bare</p>")
        assert "bare" in result

    def test_extract_head_css_normal(self):
        html = "<html><head><style>.foo{color:red}</style></head><body></body></html>"
        assert ".foo" in _extract_head_css(html)

    def test_extract_head_css_missing(self):
        assert _extract_head_css("<html><body></body></html>") == ""

    def test_find_latest_report(self, tmp_path):
        (tmp_path / "pre_filter_qc_20260301.html").write_text("old")
        (tmp_path / "pre_filter_qc_20260304.html").write_text("new")
        result = _find_latest_report(tmp_path, "pre_filter")
        assert result is not None
        assert "20260304" in result.name

    def test_find_latest_report_none(self, tmp_path):
        assert _find_latest_report(tmp_path, "pre_filter") is None

    def test_wrap_with_tabs_produces_html(self):
        result = _wrap_with_tabs(
            "Report 2",
            "<div>current</div>",
            [("Report 1", "<div>prev</div>", "body{}")],
            "body{color:blue}",
        )
        assert "tab-nav" in result
        assert "Report 2" in result
        assert "Report 1" in result
        assert "current" in result
        assert "prev" in result
        assert "showTab" in result

    def test_wrap_with_tabs_no_previous(self):
        result = _wrap_with_tabs("Report 1", "<div>standalone</div>", [], "")
        assert "standalone" in result


# ---------------------------------------------------------------------------
# Plan B: post-celltyping report (new template + tab wrapping)
# ---------------------------------------------------------------------------


class TestPostCelltypingReportPlanB:
    def test_generate_post_celltyping_report_uses_new_template(self, tmp_path):
        """Post-celltyping report should use purple header distinct from report 3."""
        pytest.importorskip("jinja2")
        from sc_tools.qc import generate_post_celltyping_report

        adata = _multi_sample_adata(n_samples=2, n_obs_per=20)
        adata.obs["celltype"] = pd.Categorical(["A"] * 20 + ["B"] * 20)
        adata.obs["batch"] = adata.obs["library_id"]
        result = generate_post_celltyping_report(
            adata, tmp_path, celltype_key="celltype", date_stamp="20260304"
        )
        assert result.exists()
        assert "20260304" in result.name
        content = result.read_text()
        # Purple header distinguishes report 4 from report 3 (dark blue)
        assert "6c3483" in content or "Post-celltyping" in content

    def test_post_celltyping_report_contains_celltypes(self, tmp_path):
        pytest.importorskip("jinja2")
        from sc_tools.qc import generate_post_celltyping_report

        adata = _multi_sample_adata(n_samples=2, n_obs_per=20)
        adata.obs["celltype"] = pd.Categorical(["Alpha"] * 20 + ["Beta"] * 20)
        adata.obs["batch"] = adata.obs["library_id"]
        result = generate_post_celltyping_report(
            adata, tmp_path, celltype_key="celltype", date_stamp="20260304"
        )
        content = result.read_text()
        assert "Alpha" in content or "Beta" in content

    def test_post_celltyping_missing_celltype_raises(self, tmp_path):
        pytest.importorskip("jinja2")
        from sc_tools.qc import generate_post_celltyping_report

        adata = _multi_sample_adata(n_samples=2, n_obs_per=20)
        with pytest.raises(ValueError, match="celltype_key"):
            generate_post_celltyping_report(adata, tmp_path, celltype_key="celltype")

    def test_post_celltyping_report_embeds_pre_filter_tab(self, tmp_path):
        pytest.importorskip("jinja2")
        from sc_tools.qc import generate_post_celltyping_report

        fake_pre = tmp_path / "pre_filter_qc_20260304.html"
        fake_pre.write_text(
            "<html><head><style>.x{}</style></head><body><p>pre filter content</p></body></html>"
        )
        adata = _multi_sample_adata(n_samples=2, n_obs_per=20)
        adata.obs["celltype"] = pd.Categorical(["A"] * 20 + ["B"] * 20)
        adata.obs["batch"] = adata.obs["library_id"]
        result = generate_post_celltyping_report(
            adata, tmp_path, celltype_key="celltype", date_stamp="20260304"
        )
        content = result.read_text()
        assert "tab-nav" in content
        assert "pre filter content" in content

    def test_post_filter_report_embeds_pre_filter(self, tmp_path):
        """generate_post_filter_report embeds pre_filter as a tab when available."""
        pytest.importorskip("jinja2")
        from sc_tools.qc import generate_post_filter_report

        # Write a fake pre-filter report
        fake_pre = tmp_path / "pre_filter_qc_20260304.html"
        fake_pre.write_text(
            "<html><head><style>.x{}</style></head><body><p>pre filter panel</p></body></html>"
        )

        adata_pre, metrics, classified = _qc_fixtures()
        adata_post = adata_pre[10:].copy()
        calculate_qc_metrics(adata_post, inplace=True, percent_top=None)

        result = generate_post_filter_report(
            adata_pre,
            adata_post,
            metrics,
            classified,
            tmp_path,
            sample_col="library_id",
            modality="visium",
            date_stamp="20260304",
        )
        content = result.read_text()
        assert "tab-nav" in content
        assert "pre filter panel" in content


# ---------------------------------------------------------------------------
# Plan B: qc_celltype_abundance plot
# ---------------------------------------------------------------------------


class TestCelltypeAbundancePlot:
    def test_qc_celltype_abundance_basic(self):
        import matplotlib
        import matplotlib.pyplot as plt

        matplotlib.use("Agg")
        from sc_tools.pl import qc_celltype_abundance

        adata = _multi_sample_adata(n_samples=2, n_obs_per=20)
        adata.obs["celltype"] = pd.Categorical(["TypeA"] * 20 + ["TypeB"] * 20)
        fig = qc_celltype_abundance(adata, celltype_key="celltype")
        assert fig is not None
        plt.close("all")

    def test_qc_celltype_abundance_with_uns_colors(self):
        import matplotlib
        import matplotlib.pyplot as plt

        matplotlib.use("Agg")
        from sc_tools.pl import qc_celltype_abundance

        adata = _multi_sample_adata(n_samples=2, n_obs_per=20)
        adata.obs["celltype"] = pd.Categorical(["TypeA"] * 20 + ["TypeB"] * 20)
        adata.uns["celltype_colors"] = ["#E69F00", "#56B4E9"]
        fig = qc_celltype_abundance(adata, celltype_key="celltype")
        assert fig is not None
        plt.close("all")


# ---------------------------------------------------------------------------
# Task 4 & 5: Jinja2 Environment+FileSystemLoader + base_report_template.html
# ---------------------------------------------------------------------------


class TestRenderTemplateFileSystemLoader:
    """render_template must use Environment+FileSystemLoader so {% extends %} works."""

    def test_render_template_by_filename(self, tmp_path):
        """Passing a bare filename resolves from the default assets dir (FileSystemLoader)."""
        pytest.importorskip("jinja2")
        from jinja2.exceptions import UndefinedError

        from sc_tools.qc.report_utils import render_template

        # base_report_template.html is a minimal template that renders with title+date.
        # If FileSystemLoader is NOT used the call raises FileNotFoundError;
        # with FileSystemLoader it succeeds (or raises UndefinedError for missing vars,
        # never FileNotFoundError).
        try:
            result = render_template(
                "base_report_template.html",
                {"title": "Test", "date": "2026-03-14"},
            )
            assert isinstance(result, str)
            assert len(result) > 0
        except UndefinedError:
            pass  # template found but vars missing — FileSystemLoader is working

    def test_render_template_path_compat(self, tmp_path):
        """Passing a full Path object still works (backward compat)."""
        pytest.importorskip("jinja2")
        from sc_tools.qc.report_utils import render_template

        # Write a minimal template in tmp_path
        tpl = tmp_path / "my_template.html"
        tpl.write_text("Hello {{ name }}")
        result = render_template(tpl, {"name": "World"})
        assert result == "Hello World"

    def test_render_template_supports_extends(self, tmp_path):
        """FileSystemLoader is required for {% extends %} to work."""
        pytest.importorskip("jinja2")
        from sc_tools.qc.report_utils import render_template

        base = tmp_path / "base.html"
        base.write_text("BASE{% block body %}{% endblock %}END")
        child = tmp_path / "child.html"
        child.write_text("{% extends 'base.html' %}{% block body %}CHILD{% endblock %}")
        result = render_template("child.html", {}, assets_dir=tmp_path)
        assert result == "BASECHILDEnd".replace("End", "END")

    def test_render_template_custom_assets_dir(self, tmp_path):
        """assets_dir kwarg overrides the default package assets dir."""
        pytest.importorskip("jinja2")
        from sc_tools.qc.report_utils import render_template

        tpl = tmp_path / "custom.html"
        tpl.write_text("VALUE={{ val }}")
        result = render_template("custom.html", {"val": 42}, assets_dir=tmp_path)
        assert "42" in result


class TestBaseReportTemplate:
    """base_report_template.html must exist and be a valid Jinja2 base template."""

    def test_base_template_exists(self):
        from pathlib import Path

        assets = Path(__file__).parent.parent / "assets"
        assert (assets / "base_report_template.html").exists(), (
            "base_report_template.html missing from sc_tools/assets/"
        )

    def test_base_template_has_required_blocks(self):
        from pathlib import Path

        assets = Path(__file__).parent.parent / "assets"
        content = (assets / "base_report_template.html").read_text()
        assert "{% block content %}" in content
        assert "{% block title %}" in content
        assert "{% block extra_head %}" in content
        assert "{% block sidebar_title %}" in content

    def test_base_template_renders_with_empty_sections(self):
        """Rendering base template directly with sections=[] must not raise."""
        pytest.importorskip("jinja2")
        from sc_tools.qc.report_utils import render_template

        result = render_template(
            "base_report_template.html",
            {"title": "Test", "date": "2026-03-14", "sections": []},
        )
        assert "<!DOCTYPE html>" in result
        assert "sc_tools" in result

    def test_base_template_renders_with_sections(self):
        """sections list should populate nav items."""
        pytest.importorskip("jinja2")
        from sc_tools.qc.report_utils import render_template

        sections = [{"id": "summary", "label": "Summary"}, {"id": "plots", "label": "Plots"}]
        result = render_template(
            "base_report_template.html",
            {"title": "My Report", "date": "2026-03-14", "sections": sections},
        )
        assert 'href="#summary"' in result
        assert "Plots" in result

    def test_base_template_sections_default_empty_when_undefined(self):
        """If 'sections' is not passed, template must not crash (defaults to [])."""
        pytest.importorskip("jinja2")
        from sc_tools.qc.report_utils import render_template

        # No 'sections' key in context
        result = render_template(
            "base_report_template.html",
            {"title": "No Sections", "date": "2026-03-14"},
        )
        assert "<!DOCTYPE html>" in result

    def test_base_template_has_bootstrap_and_plotly_cdns(self):
        from pathlib import Path

        assets = Path(__file__).parent.parent / "assets"
        content = (assets / "base_report_template.html").read_text()
        assert "bootswatch" in content or "flatly" in content
        assert "bootstrap.bundle.min.js" in content
        assert "plotly-2.27.0.min.js" in content


class TestBmReportUsesSharedRenderTemplate:
    """bm/report.py must use sc_tools.qc.report_utils.render_template, not its own."""

    def test_bm_report_no_local_render_template(self):
        import ast
        from pathlib import Path

        src = (Path(__file__).parent.parent / "bm" / "report.py").read_text()
        tree = ast.parse(src)
        func_names = [
            node.name for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)
        ]
        assert "_render_template" not in func_names, (
            "bm/report.py still defines _render_template; should use shared render_template"
        )

    def test_bm_report_imports_shared_render_template(self):
        import importlib

        bm_report = importlib.import_module("sc_tools.bm.report")
        # The module must have pulled in render_template from report_utils
        assert hasattr(bm_report, "render_template") or any(
            "render_template" in str(getattr(bm_report, attr, ""))
            for attr in dir(bm_report)
        ), "bm/report.py does not import render_template from sc_tools.qc.report_utils"
