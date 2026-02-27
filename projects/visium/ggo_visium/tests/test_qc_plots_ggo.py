"""
Project test: QC plots for ggo_visium using sc_tools.qc.

Loads adata from results/ (p2 or p35), uses an in-memory copy for QC metrics,
produces qc_2x2_grid, qc_2x4_pre_post, and qc_spatial_multipage to figures/QC/.
Does not modify or write any AnnData files; only figure files (PDF/PNG) are written.

Run from repo root: pytest projects/visium/ggo_visium/tests/test_qc_plots_ggo.py -v
Or from project dir: pytest tests/test_qc_plots_ggo.py -v (after cd projects/visium/ggo_visium)
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import pytest
import scanpy as sc

from sc_tools.qc import (
    calculate_qc_metrics,
    highly_variable_genes,
    plot_highly_variable_genes,
    plot_spatially_variable_genes,
    qc_2x2_grid,
    qc_2x4_pre_post,
    qc_scatter_counts_genes,
    qc_spatial_multipage,
    qc_violin_metrics,
    spatially_variable_genes_per_library,
)

# Project dir: parent of tests/
PROJECT_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_DIR / "results"
FIGURES_QC_RAW = PROJECT_DIR / "figures" / "QC" / "raw"
FIGURES_QC_POST = PROJECT_DIR / "figures" / "QC" / "post"

# Prefer p2 (annotated) for raw-like QC; p35 (scored) for post
ADATA_RAW_LIKE = RESULTS_DIR / "adata.annotated.p2.h5ad"
ADATA_POST = RESULTS_DIR / "adata.normalized.scored.p35.h5ad"
LIBRARY_ID_COL = "library_id"


def _ensure_qc_metrics(adata):
    """Add total_counts, n_genes_by_counts, pct_counts_mt if missing."""
    if "total_counts" in adata.obs.columns:
        return
    calculate_qc_metrics(adata, inplace=True, percent_top=(50, 100, 200, 500))


@pytest.fixture(scope="module")
def adata_raw_like():
    """Load p2 if present; skip test if missing. Uses a copy so no .h5ad is modified or written."""
    if not ADATA_RAW_LIKE.exists():
        pytest.skip(f"Missing {ADATA_RAW_LIKE}; run Phase 1–2 first.")
    adata = sc.read_h5ad(ADATA_RAW_LIKE)
    adata = adata.copy()
    _ensure_qc_metrics(adata)
    return adata


@pytest.fixture(scope="module")
def adata_post():
    """Load p35 if present; skip if missing. Uses a copy so no .h5ad is modified or written."""
    if not ADATA_POST.exists():
        pytest.skip(f"Missing {ADATA_POST}; run Phase 3.5b first.")
    adata = sc.read_h5ad(ADATA_POST)
    adata = adata.copy()
    _ensure_qc_metrics(adata)
    return adata


def test_qc_2x2_raw_ggo_visium(adata_raw_like):
    """Produce 2x2 QC grid for raw-like (p2) and save to figures/QC/raw."""
    FIGURES_QC_RAW.mkdir(parents=True, exist_ok=True)
    fig = qc_2x2_grid(
        adata_raw_like,
        output_dir=FIGURES_QC_RAW,
        basename="qc_2x2_raw",
        dpi=300,
    )
    assert fig is not None
    assert (FIGURES_QC_RAW / "qc_2x2_raw.pdf").exists()
    assert (FIGURES_QC_RAW / "qc_2x2_raw.png").exists()
    fig.clear()
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_qc_2x2_post_ggo_visium(adata_post):
    """Produce 2x2 QC grid for post (p35) and save to figures/QC/post."""
    FIGURES_QC_POST.mkdir(parents=True, exist_ok=True)
    fig = qc_2x2_grid(
        adata_post,
        output_dir=FIGURES_QC_POST,
        basename="qc_2x2_post",
        dpi=300,
    )
    assert fig is not None
    assert (FIGURES_QC_POST / "qc_2x2_post.pdf").exists()
    assert (FIGURES_QC_POST / "qc_2x2_post.png").exists()
    fig.clear()
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_qc_2x4_pre_post_ggo_visium(adata_raw_like, adata_post):
    """Produce 2x4 pre vs post QC (pre 2x2 left, post 2x2 right) and save to figures/QC/post."""
    FIGURES_QC_POST.mkdir(parents=True, exist_ok=True)
    fig = qc_2x4_pre_post(
        adata_raw_like,
        adata_post,
        output_dir=FIGURES_QC_POST,
        basename="qc_2x4_pre_post",
        dpi=300,
    )
    assert fig is not None
    assert (FIGURES_QC_POST / "qc_2x4_pre_post.pdf").exists()
    assert (FIGURES_QC_POST / "qc_2x4_pre_post.png").exists()
    fig.clear()
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_qc_spatial_multipage_raw_ggo_visium(adata_raw_like):
    """Produce multipage spatial QC (total_count, log1p, %mt) for raw-like; save to figures/QC/raw."""
    if LIBRARY_ID_COL not in adata_raw_like.obs.columns:
        pytest.skip(f"obs lacks {LIBRARY_ID_COL}; required for multipage.")
    if "spatial" not in adata_raw_like.obsm:
        pytest.skip("obsm['spatial'] required for multipage.")
    spatial_uns = adata_raw_like.uns.get("spatial", {})
    libs = adata_raw_like.obs[LIBRARY_ID_COL].dropna().unique()
    if not spatial_uns or not all(str(l) in spatial_uns for l in libs):
        pytest.skip("uns['spatial'] must contain keys for each library_id for multipage.")
    FIGURES_QC_RAW.mkdir(parents=True, exist_ok=True)
    out_pdf = FIGURES_QC_RAW / "qc_spatial_multipage_raw.pdf"
    qc_spatial_multipage(
        adata_raw_like,
        LIBRARY_ID_COL,
        str(out_pdf),
        total_counts_col="total_counts",
        pct_mt_col="pct_counts_mt",
        common_scale=True,
    )
    assert out_pdf.exists()


def test_qc_spatial_multipage_post_ggo_visium(adata_post):
    """Produce multipage spatial QC for post (p35); save to figures/QC/post."""
    if LIBRARY_ID_COL not in adata_post.obs.columns:
        pytest.skip(f"obs lacks {LIBRARY_ID_COL}; required for multipage.")
    if "spatial" not in adata_post.obsm:
        pytest.skip("obsm['spatial'] required for multipage.")
    spatial_uns = adata_post.uns.get("spatial", {})
    libs = adata_post.obs[LIBRARY_ID_COL].dropna().unique()
    if not spatial_uns or not all(str(l) in spatial_uns for l in libs):
        pytest.skip("uns['spatial'] must contain keys for each library_id for multipage.")
    FIGURES_QC_POST.mkdir(parents=True, exist_ok=True)
    out_pdf = FIGURES_QC_POST / "qc_spatial_multipage_post.pdf"
    qc_spatial_multipage(
        adata_post,
        LIBRARY_ID_COL,
        str(out_pdf),
        total_counts_col="total_counts",
        pct_mt_col="pct_counts_mt",
        common_scale=True,
    )
    assert out_pdf.exists()


def test_qc_violin_raw_ggo_visium(adata_raw_like):
    """Produce violin QC metrics for raw-like; save to figures/QC/raw."""
    FIGURES_QC_RAW.mkdir(parents=True, exist_ok=True)
    fig = qc_violin_metrics(
        adata_raw_like,
        output_dir=FIGURES_QC_RAW,
        basename="qc_violin_raw",
        dpi=300,
    )
    assert fig is not None
    assert (FIGURES_QC_RAW / "qc_violin_raw.pdf").exists()
    fig.clear()
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_qc_violin_post_ggo_visium(adata_post):
    """Produce violin QC metrics for post; save to figures/QC/post."""
    FIGURES_QC_POST.mkdir(parents=True, exist_ok=True)
    fig = qc_violin_metrics(
        adata_post,
        output_dir=FIGURES_QC_POST,
        basename="qc_violin_post",
        dpi=300,
    )
    assert fig is not None
    fig.clear()
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_qc_scatter_raw_ggo_visium(adata_raw_like):
    """Produce scatter counts vs genes for raw-like; save to figures/QC/raw."""
    FIGURES_QC_RAW.mkdir(parents=True, exist_ok=True)
    fig = qc_scatter_counts_genes(
        adata_raw_like,
        output_dir=FIGURES_QC_RAW,
        basename="qc_scatter_raw",
        dpi=300,
    )
    assert fig is not None
    assert (FIGURES_QC_RAW / "qc_scatter_raw.pdf").exists()
    fig.clear()
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_qc_scatter_post_ggo_visium(adata_post):
    """Produce scatter counts vs genes for post; save to figures/QC/post."""
    FIGURES_QC_POST.mkdir(parents=True, exist_ok=True)
    fig = qc_scatter_counts_genes(
        adata_post,
        output_dir=FIGURES_QC_POST,
        basename="qc_scatter_post",
        dpi=300,
    )
    assert fig is not None
    fig.clear()
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_plot_highly_variable_genes_ggo_visium(adata_post):
    """Produce HVG plot for post; run highly_variable_genes on copy if var has no HVG."""
    FIGURES_QC_POST.mkdir(parents=True, exist_ok=True)
    adata = adata_post
    if "highly_variable" not in adata.var.columns:
        adata = adata.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000, inplace=True)
    fig = plot_highly_variable_genes(
        adata,
        output_dir=FIGURES_QC_POST,
        basename="hvg_post",
        dpi=300,
    )
    assert fig is not None
    assert (FIGURES_QC_POST / "hvg_post.pdf").exists()
    fig.clear()
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_plot_spatially_variable_genes_ggo_visium(adata_post):
    """Produce SVG plot when library_id present; skip if no library_id (no adata write)."""
    if LIBRARY_ID_COL not in adata_post.obs.columns:
        pytest.skip(f"obs lacks {LIBRARY_ID_COL}; skip SVG.")
    if "spatial" not in adata_post.obsm:
        pytest.skip("obsm['spatial'] required for SVG.")
    per_lib = spatially_variable_genes_per_library(
        adata_post,
        library_id_col=LIBRARY_ID_COL,
        genes=adata_post.var_names[: min(100, adata_post.n_vars)].tolist(),
        n_top_genes=20,
    )
    if per_lib is None:
        pytest.skip("spatially_variable_genes_per_library returned None (e.g. no squidpy or too few obs).")
    FIGURES_QC_POST.mkdir(parents=True, exist_ok=True)
    fig = plot_spatially_variable_genes(
        adata_post,
        output_dir=FIGURES_QC_POST,
        basename="svg_post",
        dpi=300,
    )
    assert fig is not None
    assert (FIGURES_QC_POST / "svg_post.pdf").exists()
    fig.clear()
    import matplotlib.pyplot as plt
    plt.close(fig)
