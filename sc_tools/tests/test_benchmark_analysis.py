"""Tests for sc_tools.bm.analysis — cross-dataset statistics."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from sc_tools.bm.analysis import (
    compute_generalization_matrix,
    identify_failure_cases,
    rank_methods_per_tissue,
    statistical_comparison,
)


def _make_results_df(n_rois: int = 20) -> pd.DataFrame:
    """Create synthetic benchmark results."""
    rng = np.random.RandomState(42)
    methods = ["cellpose_cyto2", "stardist", "deepcell"]
    datasets = ["ggo-imc", "aapc"]
    tissues = ["lung", "colon"]

    rows = []
    for i in range(n_rois):
        for method in methods:
            rows.append(
                {
                    "roi_id": f"roi_{i}",
                    "method": method,
                    "dataset": datasets[i % 2],
                    "tissue": tissues[i % 2],
                    "strategy": 1,
                    "n_cells": rng.randint(50, 500),
                    "boundary_regularity": rng.uniform(0.5, 0.95),
                    "median_area": rng.uniform(50, 200),
                    "runtime_s": rng.uniform(1, 10),
                }
            )
    return pd.DataFrame(rows)


class TestGeneralizationMatrix:
    def test_basic(self):
        df = _make_results_df()
        matrix = compute_generalization_matrix(df, metric="boundary_regularity")
        assert isinstance(matrix, pd.DataFrame)
        assert "overall_mean" in matrix.columns
        assert len(matrix) == 3  # 3 methods

    def test_invalid_metric(self):
        df = _make_results_df()
        with pytest.raises(ValueError):
            compute_generalization_matrix(df, metric="nonexistent")


class TestRankMethods:
    def test_ranking(self):
        df = _make_results_df()
        ranked = rank_methods_per_tissue(df, metric="boundary_regularity")
        assert "rank" in ranked.columns
        assert "tissue" in ranked.columns
        # Each tissue should have all methods
        for tissue in df["tissue"].unique():
            tissue_ranked = ranked[ranked["tissue"] == tissue]
            assert len(tissue_ranked) == df["method"].nunique()


class TestStatisticalComparison:
    def test_pairwise_wilcoxon(self):
        df = _make_results_df(n_rois=30)
        result = statistical_comparison(df, metric="boundary_regularity")
        assert isinstance(result, pd.DataFrame)
        if len(result) > 0:
            assert "p_value" in result.columns
            assert "p_adjusted" in result.columns
            assert "significant" in result.columns
            # BH correction: adjusted >= raw
            assert all(result["p_adjusted"] >= result["p_value"] - 1e-10)

    def test_no_roi_id(self):
        df = _make_results_df()
        df = df.drop(columns=["roi_id"])
        result = statistical_comparison(df, metric="boundary_regularity")
        assert len(result) == 0


class TestFailureCases:
    def test_identify_worst(self):
        df = _make_results_df()
        worst = identify_failure_cases(df, metric="boundary_regularity", n_worst=3)
        assert isinstance(worst, pd.DataFrame)
        # Each method should have at most n_worst entries
        for method in df["method"].unique():
            method_worst = worst[worst["method"] == method]
            assert len(method_worst) <= 3
