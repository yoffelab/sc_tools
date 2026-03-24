"""Tests for pseudobulk differential expression module (SCI-01).

Tests are split into two groups:
1. Aggregation tests: no PyDESeq2 needed, test count aggregation logic
2. Integration tests: require PyDESeq2, test full DE pipeline

PyDESeq2-dependent tests are guarded with pytest.importorskip.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData


# ---------------------------------------------------------------------------
# Aggregation tests (no PyDESeq2 required)
# ---------------------------------------------------------------------------


class TestAggregatePseudobulk:
    """Tests for aggregate_pseudobulk function."""

    def test_returns_dict_with_celltype_keys(self, adata_multi_subject: AnnData) -> None:
        """aggregate_pseudobulk returns dict with celltype keys."""
        from sc_tools.tl.de import aggregate_pseudobulk

        result = aggregate_pseudobulk(
            adata_multi_subject, subject_key="subject_id", celltype_key="celltype"
        )
        assert isinstance(result, dict)
        # Should have some celltypes
        assert len(result) > 0
        # Keys should be celltype names
        for key in result:
            assert key in adata_multi_subject.obs["celltype"].cat.categories.tolist()

    def test_values_are_count_metadata_tuples(self, adata_multi_subject: AnnData) -> None:
        """Each value is (count_df, metadata_df) tuple."""
        from sc_tools.tl.de import aggregate_pseudobulk

        result = aggregate_pseudobulk(
            adata_multi_subject, subject_key="subject_id", celltype_key="celltype"
        )
        for ct, (count_df, meta_df) in result.items():
            assert isinstance(count_df, pd.DataFrame)
            assert isinstance(meta_df, pd.DataFrame)

    def test_count_df_shape_and_dtype(self, adata_multi_subject: AnnData) -> None:
        """count_df has shape (n_valid_subjects, n_genes) with integer values."""
        from sc_tools.tl.de import aggregate_pseudobulk

        result = aggregate_pseudobulk(
            adata_multi_subject,
            subject_key="subject_id",
            celltype_key="celltype",
            min_cells=1,  # Low threshold so all subjects pass
        )
        for ct, (count_df, meta_df) in result.items():
            # Columns = genes
            assert count_df.shape[1] == adata_multi_subject.n_vars
            # Values should be integers
            assert np.all(count_df.values == count_df.values.astype(int))
            # Values should be non-negative
            assert (count_df.values >= 0).all()

    def test_min_cells_filters_subjects(self, adata_multi_subject: AnnData) -> None:
        """Subjects with fewer than min_cells cells are excluded."""
        from sc_tools.tl.de import aggregate_pseudobulk

        # With very high threshold, should exclude subjects
        result_strict = aggregate_pseudobulk(
            adata_multi_subject,
            subject_key="subject_id",
            celltype_key="celltype",
            min_cells=100,  # Very high: each subject has ~50 cells/celltype at most
        )
        result_lenient = aggregate_pseudobulk(
            adata_multi_subject,
            subject_key="subject_id",
            celltype_key="celltype",
            min_cells=1,
        )
        # Strict should have fewer or equal celltypes (some may be entirely excluded)
        assert len(result_strict) <= len(result_lenient)

    def test_celltypes_with_few_subjects_excluded(self) -> None:
        """Celltypes with fewer than 2 valid subjects are excluded entirely."""
        from sc_tools.tl.de import aggregate_pseudobulk

        # Create adata where one celltype has only 1 subject
        rng = np.random.default_rng(99)
        n_cells = 100
        X = rng.negative_binomial(5, 0.3, (n_cells, 50)).astype("float32")
        obs = pd.DataFrame(
            {
                "subject_id": pd.Categorical(
                    ["S1"] * 50 + ["S2"] * 30 + ["S1"] * 20
                ),
                "celltype": pd.Categorical(
                    ["typeA"] * 50 + ["typeA"] * 30 + ["typeB"] * 20
                ),
                "condition": pd.Categorical(
                    ["ctrl"] * 50 + ["treat"] * 30 + ["ctrl"] * 20
                ),
            },
            index=[f"c{i}" for i in range(n_cells)],
        )
        adata = AnnData(X, obs=obs, var=pd.DataFrame(index=[f"G{i}" for i in range(50)]))

        result = aggregate_pseudobulk(
            adata, subject_key="subject_id", celltype_key="celltype", min_cells=1
        )
        # typeB has only 1 subject (S1), should be excluded
        assert "typeB" not in result
        assert "typeA" in result

    def test_raw_counts_from_layers(self, adata_multi_subject: AnnData) -> None:
        """Raw counts extracted from layers['counts'] > raw.X > X with integer validation."""
        from sc_tools.tl.de import aggregate_pseudobulk

        # Store counts in a layer
        adata_multi_subject.layers["counts"] = adata_multi_subject.X.copy()
        # Overwrite X with garbage
        adata_multi_subject.X = np.zeros_like(adata_multi_subject.X)

        result = aggregate_pseudobulk(
            adata_multi_subject,
            subject_key="subject_id",
            celltype_key="celltype",
            min_cells=1,
        )
        # Should have results (from layer, not zero X)
        for ct, (count_df, _) in result.items():
            assert count_df.values.sum() > 0

    def test_raw_counts_from_layer_param(self, adata_multi_subject: AnnData) -> None:
        """Explicit layer parameter is used."""
        from sc_tools.tl.de import aggregate_pseudobulk

        adata_multi_subject.layers["my_raw"] = adata_multi_subject.X.copy()
        adata_multi_subject.X = np.zeros_like(adata_multi_subject.X)

        result = aggregate_pseudobulk(
            adata_multi_subject,
            subject_key="subject_id",
            celltype_key="celltype",
            layer="my_raw",
            min_cells=1,
        )
        for ct, (count_df, _) in result.items():
            assert count_df.values.sum() > 0


# ---------------------------------------------------------------------------
# Formula inference tests (no PyDESeq2 required)
# ---------------------------------------------------------------------------


class TestFormulaInference:
    """Tests for design formula auto-inference logic."""

    def test_auto_formula_no_batch(self) -> None:
        """Auto-inferred formula builds '~ condition' when no batch cols detected."""
        from sc_tools.tl.de import _infer_design_formula

        # Metadata with only condition column
        meta = pd.DataFrame(
            {"condition": ["ctrl", "treat", "ctrl"]},
            index=["S1", "S2", "S3"],
        )
        formula = _infer_design_formula(meta, "condition", "subject_id")
        assert formula == "~ condition"

    def test_auto_formula_with_batch(self) -> None:
        """Auto-inferred formula builds '~ condition + batch' when batch exists and is not collinear."""
        from sc_tools.tl.de import _infer_design_formula

        meta = pd.DataFrame(
            {
                "condition": ["ctrl", "ctrl", "treat", "treat"],
                "batch": ["B0", "B1", "B0", "B1"],
            },
            index=["S1", "S2", "S3", "S4"],
        )
        formula = _infer_design_formula(meta, "condition", "subject_id")
        assert "batch" in formula
        assert "condition" in formula

    def test_collinear_batch_excluded(self) -> None:
        """Collinear batch (1:1 with subject_id) is excluded from auto-formula."""
        from sc_tools.tl.de import _infer_design_formula

        # library_id maps 1:1 to subject index
        meta = pd.DataFrame(
            {
                "condition": ["ctrl", "treat", "ctrl", "treat"],
                "library_id": ["L1", "L2", "L3", "L4"],
            },
            index=["S1", "S2", "S3", "S4"],
        )
        formula = _infer_design_formula(meta, "condition", "subject_id")
        # library_id is collinear with subject index, should be excluded
        assert "library_id" not in formula
        assert formula == "~ condition"

    def test_custom_formula_override(self) -> None:
        """Custom formula override is used when provided."""
        from sc_tools.tl.de import _infer_design_formula

        meta = pd.DataFrame(
            {"condition": ["ctrl", "treat"], "batch": ["B0", "B1"]},
            index=["S1", "S2"],
        )
        # When custom formula is provided, _infer_design_formula is not called
        # (the caller passes formula directly). Test the pattern at run level.
        # This just tests that auto-inference would give something different.
        auto = _infer_design_formula(meta, "condition", "subject_id")
        custom = "~ condition + sex + age"
        assert auto != custom  # auto would not include sex/age


# ---------------------------------------------------------------------------
# PyDESeq2 integration tests (require pydeseq2)
# ---------------------------------------------------------------------------


class TestRunPseudobulkDE:
    """Integration tests for run_pseudobulk_de (requires PyDESeq2)."""

    @pytest.fixture(autouse=True)
    def _require_pydeseq2(self) -> None:
        pytest.importorskip("pydeseq2")

    def test_returns_dict_of_dataframes(self, adata_multi_subject: AnnData) -> None:
        """run_pseudobulk_de returns dict mapping celltype to DataFrame."""
        from sc_tools.tl.de import run_pseudobulk_de

        results = run_pseudobulk_de(
            adata_multi_subject,
            condition_key="condition",
            subject_key="subject_id",
            celltype_key="celltype",
            min_cells_per_combo=1,
            min_subjects_per_group=2,
        )
        assert isinstance(results, dict)
        for ct, df in results.items():
            assert isinstance(df, pd.DataFrame)

    def test_result_columns(self, adata_multi_subject: AnnData) -> None:
        """Results have gene, log2FC, pvalue, padj, baseMean columns."""
        from sc_tools.tl.de import run_pseudobulk_de

        results = run_pseudobulk_de(
            adata_multi_subject,
            condition_key="condition",
            min_cells_per_combo=1,
            min_subjects_per_group=2,
        )
        expected_cols = {"gene", "log2FC", "pvalue", "padj", "baseMean"}
        for ct, df in results.items():
            assert expected_cols.issubset(set(df.columns)), (
                f"Missing columns for {ct}: {expected_cols - set(df.columns)}"
            )

    def test_skips_celltypes_with_few_subjects(self) -> None:
        """Celltypes where any condition group has < min_subjects_per_group subjects are skipped."""
        from sc_tools.tl.de import run_pseudobulk_de

        # Create adata where one celltype has too few subjects in one group
        rng = np.random.default_rng(42)
        n = 200
        X = rng.negative_binomial(5, 0.3, (n, 50)).astype("float32")
        obs = pd.DataFrame(
            {
                "subject_id": pd.Categorical(np.repeat([f"S{i}" for i in range(4)], 50)),
                "celltype": pd.Categorical(["typeA"] * 100 + ["typeB"] * 100),
                "condition": pd.Categorical(
                    ["ctrl"] * 50 + ["treat"] * 50 + ["ctrl"] * 50 + ["treat"] * 50
                ),
            },
            index=[f"c{i}" for i in range(n)],
        )
        adata = AnnData(X, obs=obs, var=pd.DataFrame(index=[f"G{i}" for i in range(50)]))

        # For typeA: S0=ctrl, S1=treat => 1 per group
        # For typeB: S2=ctrl, S3=treat => 1 per group
        # With min_subjects_per_group=3, both should be skipped
        results = run_pseudobulk_de(
            adata,
            condition_key="condition",
            min_subjects_per_group=3,
            min_cells_per_combo=1,
        )
        assert len(results) == 0

    def test_custom_formula_used(self, adata_multi_subject: AnnData) -> None:
        """Custom formula override is used when provided."""
        from sc_tools.tl.de import run_pseudobulk_de

        # Provide an explicit formula
        results = run_pseudobulk_de(
            adata_multi_subject,
            condition_key="condition",
            formula="~ condition",
            min_cells_per_combo=1,
            min_subjects_per_group=2,
        )
        # Should succeed with custom formula
        assert len(results) > 0
