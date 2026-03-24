"""Tests for marker validation compute and report integration (SCI-02)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def adata_with_markers() -> AnnData:
    """AnnData with 3 celltypes and known marker expression.

    - T_cell: high expression of CD3D, CD3E, CD4 (above threshold)
    - B_cell: high expression of CD19, MS4A1, CD79A (above threshold)
    - Macrophage: very low expression of all its markers (below 0.1 threshold)
    """
    rng = np.random.default_rng(42)
    n_cells = 300
    n_genes = 20

    # Gene names including markers
    gene_names = [
        "CD3D", "CD3E", "CD4", "CD8A",  # T cell markers
        "CD19", "MS4A1", "CD79A",  # B cell markers
        "CD68", "CSF1R", "CD14",  # Macrophage markers
        "ACTB", "GAPDH", "MALAT1", "RPL13",  # housekeeping
        "GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_E", "GENE_F",
    ]

    X = rng.random((n_cells, n_genes)).astype(np.float32) * 0.01  # very low baseline

    celltypes = (["T_cell"] * 100) + (["B_cell"] * 100) + (["Macrophage"] * 100)

    # T_cell markers: high expression in T_cell subset (rows 0-99)
    X[0:100, 0] = rng.random(100) * 2 + 1.0  # CD3D
    X[0:100, 1] = rng.random(100) * 2 + 1.0  # CD3E
    X[0:100, 2] = rng.random(100) * 1.5 + 0.5  # CD4

    # B_cell markers: high expression in B_cell subset (rows 100-199)
    X[100:200, 4] = rng.random(100) * 2 + 1.0  # CD19
    X[100:200, 5] = rng.random(100) * 2 + 1.0  # MS4A1
    X[100:200, 6] = rng.random(100) * 1.5 + 0.5  # CD79A

    # Macrophage markers: keep very low (below 0.1 threshold)
    # Already set to ~0.01 from baseline

    adata = AnnData(
        X=X,
        obs=pd.DataFrame(
            {"celltype": pd.Categorical(celltypes)},
            index=[f"cell_{i}" for i in range(n_cells)],
        ),
        var=pd.DataFrame(index=gene_names),
    )
    return adata


@pytest.fixture()
def marker_genes() -> dict[str, list[str]]:
    """Marker genes dict mapping celltype to marker list."""
    return {
        "T_cell": ["CD3D", "CD3E", "CD4"],
        "B_cell": ["CD19", "MS4A1", "CD79A"],
        "Macrophage": ["CD68", "CSF1R", "CD14"],
    }


# ---------------------------------------------------------------------------
# Unit tests: compute_marker_validation
# ---------------------------------------------------------------------------


class TestComputeMarkerValidation:
    """Tests for compute_marker_validation function."""

    def test_returns_dataframe_with_correct_columns(
        self, adata_with_markers, marker_genes
    ):
        from sc_tools.qc.marker_validation import compute_marker_validation

        df, summary = compute_marker_validation(
            adata_with_markers, "celltype", marker_genes
        )
        assert isinstance(df, pd.DataFrame)
        assert set(df.columns) >= {"celltype", "marker_gene", "mean_expr", "flagged"}

    def test_low_expression_type_flagged(self, adata_with_markers, marker_genes):
        """Macrophage has all markers below threshold -> flagged=True."""
        from sc_tools.qc.marker_validation import compute_marker_validation

        df, _ = compute_marker_validation(
            adata_with_markers, "celltype", marker_genes, threshold=0.1
        )
        mac_rows = df[df["celltype"] == "Macrophage"]
        assert len(mac_rows) == 3
        assert mac_rows["flagged"].all(), "Macrophage should be flagged (all markers below threshold)"

    def test_high_expression_type_not_flagged(self, adata_with_markers, marker_genes):
        """T_cell has markers above threshold -> flagged=False."""
        from sc_tools.qc.marker_validation import compute_marker_validation

        df, _ = compute_marker_validation(
            adata_with_markers, "celltype", marker_genes, threshold=0.1
        )
        tcell_rows = df[df["celltype"] == "T_cell"]
        assert not tcell_rows["flagged"].any(), "T_cell should not be flagged"

    def test_missing_genes_have_nan(self, adata_with_markers):
        """Genes not in adata.var_names should have NaN mean_expr."""
        from sc_tools.qc.marker_validation import compute_marker_validation

        marker_genes_with_missing = {
            "T_cell": ["CD3D", "NONEXISTENT_GENE"],
        }
        df, _ = compute_marker_validation(
            adata_with_markers, "celltype", marker_genes_with_missing
        )
        missing_row = df[df["marker_gene"] == "NONEXISTENT_GENE"]
        assert len(missing_row) == 1
        assert pd.isna(missing_row["mean_expr"].iloc[0])

    def test_empty_marker_genes_returns_empty_df(self, adata_with_markers):
        """Empty marker_genes dict returns empty DataFrame."""
        from sc_tools.qc.marker_validation import compute_marker_validation

        df, summary = compute_marker_validation(adata_with_markers, "celltype", {})
        assert len(df) == 0
        assert summary["n_types_tested"] == 0

    def test_summary_stats(self, adata_with_markers, marker_genes):
        """summary_stats returns dict with n_types_tested, n_flagged, total_cells."""
        from sc_tools.qc.marker_validation import compute_marker_validation

        _, summary = compute_marker_validation(
            adata_with_markers, "celltype", marker_genes, threshold=0.1
        )
        assert isinstance(summary, dict)
        assert summary["n_types_tested"] == 3
        assert summary["n_flagged"] == 1  # only Macrophage
        assert summary["total_cells"] == 300


# ---------------------------------------------------------------------------
# Unit tests: render_marker_dotplot
# ---------------------------------------------------------------------------


class TestRenderMarkerDotplot:
    """Tests for render_marker_dotplot function."""

    def test_returns_base64_string(self, adata_with_markers, marker_genes):
        from sc_tools.qc.marker_validation import render_marker_dotplot

        result = render_marker_dotplot(
            adata_with_markers, "celltype", marker_genes
        )
        assert isinstance(result, str)
        # base64 encoded PNG should be non-empty
        if result:  # may be empty if scanpy not available
            assert len(result) > 100  # non-trivial base64 string
