"""
Unit tests for sc_tools.tl.gsea.

Tests run_ora and run_gsea_pseudobulk with synthetic AnnData fixtures.
run_gsea_pseudobulk tests are skipped if gseapy is not installed.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from sc_tools.tl.gsea import run_gsea_pseudobulk, run_ora

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _make_adata(n_obs: int = 80, n_vars: int = 60, seed: int = 42):
    """Build synthetic adata with obs['group'] and obs['celltype']."""
    rng = np.random.default_rng(seed)
    # Two groups with different expression profiles
    half = n_obs // 2
    X_a = rng.negative_binomial(10, 0.3, (half, n_vars)).astype(np.float32)
    X_b = rng.negative_binomial(2, 0.5, (half, n_vars)).astype(np.float32)
    X = np.vstack([X_a, X_b])
    var_names = [f"GENE{i}" for i in range(1, n_vars + 1)]
    obs = pd.DataFrame(
        {
            "group": ["A"] * half + ["B"] * half,
            "celltype": ["TypeX"] * (half // 2)
            + ["TypeY"] * (half // 2)
            + ["TypeX"] * (half // 2)
            + ["TypeY"] * (half // 2),
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    adata = sc.AnnData(X, obs=obs, var=pd.DataFrame(index=var_names))
    return adata


# ---------------------------------------------------------------------------
# run_ora tests
# ---------------------------------------------------------------------------


def test_run_ora_returns_dataframe():
    adata = _make_adata()
    gene_sets = {"CatA": {"SetA": ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"]}}
    result = run_ora(adata, groupby="group", gene_set_dict=gene_sets)
    assert isinstance(result, pd.DataFrame)


def test_run_ora_columns():
    adata = _make_adata()
    gene_sets = {"CatA": {"SetA": [f"GENE{i}" for i in range(1, 10)]}}
    result = run_ora(adata, groupby="group", gene_set_dict=gene_sets)
    expected_cols = {
        "group",
        "category",
        "gene_set",
        "n_overlap",
        "n_set",
        "n_background",
        "n_group_genes",
        "odds_ratio",
        "p_val",
        "p_adj",
    }
    assert expected_cols.issubset(set(result.columns))


def test_run_ora_groups_present():
    adata = _make_adata()
    gene_sets = {"Cat": {"Sig": [f"GENE{i}" for i in range(1, 10)]}}
    result = run_ora(adata, groupby="group", gene_set_dict=gene_sets)
    assert set(result["group"].unique()) == {"A", "B"}


def test_run_ora_p_adj_range():
    adata = _make_adata()
    gene_sets = {"Cat": {"Sig": [f"GENE{i}" for i in range(1, 10)]}}
    result = run_ora(adata, groupby="group", gene_set_dict=gene_sets)
    assert (result["p_adj"] >= 0).all() and (result["p_adj"] <= 1).all()


def test_run_ora_flat_gene_set_dict():
    """Accept flat {name: [genes]} format."""
    adata = _make_adata()
    gene_sets = {"SetA": ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"]}
    result = run_ora(adata, groupby="group", gene_set_dict=gene_sets)
    assert len(result) >= 1


def test_run_ora_multiple_categories():
    adata = _make_adata()
    gene_sets = {
        "CatA": {"Sig1": [f"GENE{i}" for i in range(1, 8)]},
        "CatB": {"Sig2": [f"GENE{i}" for i in range(10, 20)]},
    }
    result = run_ora(adata, groupby="group", gene_set_dict=gene_sets)
    assert set(result["category"].unique()).issubset({"CatA", "CatB"})


def test_run_ora_below_min_genes_excluded():
    """Sets with fewer than min_genes genes should be excluded."""
    adata = _make_adata()
    gene_sets = {"Cat": {"TinySet": ["GENE1", "GENE2"]}}  # 2 < min_genes=3
    result = run_ora(adata, groupby="group", gene_set_dict=gene_sets, min_genes=3)
    assert len(result) == 0


def test_run_ora_custom_background():
    adata = _make_adata()
    bg = [f"GENE{i}" for i in range(1, 30)]
    gene_sets = {"Cat": {"Sig": [f"GENE{i}" for i in range(1, 10)]}}
    result = run_ora(adata, groupby="group", gene_set_dict=gene_sets, background=bg)
    assert result.iloc[0]["n_background"] == len({g.upper() for g in bg})


def test_run_ora_invalid_groupby():
    adata = _make_adata()
    with pytest.raises(ValueError, match="groupby="):
        run_ora(adata, groupby="nonexistent", gene_set_dict={"Cat": {"S": ["G1"]}})


def test_run_ora_empty_result_when_no_valid_sets():
    adata = _make_adata()
    # All sets below min_genes
    gene_sets = {"Cat": {"S": ["GENE1"]}}
    result = run_ora(adata, groupby="group", gene_set_dict=gene_sets, min_genes=5)
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 0


def test_run_ora_sorted_by_p_adj():
    adata = _make_adata(n_obs=100, n_vars=80)
    gene_sets = {"Cat": {f"Sig{i}": [f"GENE{j}" for j in range(i, i + 8)] for i in range(1, 10)}}
    result = run_ora(adata, groupby="group", gene_set_dict=gene_sets)
    for grp, sub in result.groupby("group"):
        assert (sub["p_adj"].diff().dropna() >= 0).all(), f"Not sorted for group {grp}"


# ---------------------------------------------------------------------------
# run_gsea_pseudobulk tests
# ---------------------------------------------------------------------------


def test_run_gsea_pseudobulk_import_error():
    """run_gsea_pseudobulk raises ImportError with install hint if gseapy missing."""
    try:
        import gseapy  # noqa: F401

        pytest.skip("gseapy is installed; cannot test missing-import path")
    except ImportError:
        pass

    adata = _make_adata()
    with pytest.raises(ImportError, match="gseapy"):
        run_gsea_pseudobulk(adata, groupby="group", gene_set_dict={"Cat": {"Sig": ["GENE1"]}})


def test_run_gsea_pseudobulk_returns_dataframe():
    pytest.importorskip("gseapy")
    adata = _make_adata(n_obs=60, n_vars=80)
    gene_sets = {"Cat": {"Sig": [f"GENE{i}" for i in range(1, 20)]}}
    result = run_gsea_pseudobulk(
        adata, groupby="group", gene_set_dict=gene_sets, n_permutations=100
    )
    assert isinstance(result, pd.DataFrame)


def test_run_gsea_pseudobulk_columns():
    pytest.importorskip("gseapy")
    adata = _make_adata(n_obs=60, n_vars=80)
    gene_sets = {"Cat": {"Sig": [f"GENE{i}" for i in range(1, 20)]}}
    result = run_gsea_pseudobulk(
        adata, groupby="group", gene_set_dict=gene_sets, n_permutations=100
    )
    expected = {"group", "category", "gene_set", "NES", "p_val", "p_adj", "lead_edge_genes"}
    assert expected.issubset(set(result.columns))


def test_run_gsea_pseudobulk_invalid_groupby():
    pytest.importorskip("gseapy")
    adata = _make_adata()
    with pytest.raises(ValueError, match="groupby="):
        run_gsea_pseudobulk(adata, groupby="bad_col", gene_set_dict={"Cat": {"Sig": ["GENE1"]}})


def test_run_gsea_pseudobulk_ranking_stats():
    pytest.importorskip("gseapy")
    adata = _make_adata(n_obs=60, n_vars=80)
    gene_sets = {"Cat": {"Sig": [f"GENE{i}" for i in range(1, 20)]}}
    for stat in ("logfc", "zscore", "tstat"):
        result = run_gsea_pseudobulk(
            adata, groupby="group", gene_set_dict=gene_sets, ranking_stat=stat, n_permutations=50
        )
        assert isinstance(result, pd.DataFrame), f"Failed for ranking_stat={stat}"


def test_run_gsea_pseudobulk_invalid_ranking_stat():
    pytest.importorskip("gseapy")
    adata = _make_adata()
    with pytest.raises(ValueError, match="ranking_stat="):
        run_gsea_pseudobulk(
            adata,
            groupby="group",
            gene_set_dict={"Cat": {"Sig": [f"GENE{i}" for i in range(20)]}},
            ranking_stat="bad",
        )
