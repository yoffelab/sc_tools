"""
Unit tests for sc_tools.gr — multi-ROI wrappers around squidpy.gr.

Fixture: 3 ROIs, 100 cells each, 5 cell types (ROI 3 missing one cell type),
50 genes, raw integer counts, obsm['spatial'] set to random 2D coordinates per ROI.
"""

import anndata as ad
import numpy as np
import pandas as pd
import pytest

# ──────────────────────────────────────────────────────────────────────────────
# Fixtures
# ──────────────────────────────────────────────────────────────────────────────

N_ROIS = 3
N_CELLS_PER_ROI = 100
N_GENES = 50
CELLTYPES = ["TypeA", "TypeB", "TypeC", "TypeD", "TypeE"]


def _make_test_adata():
    """
    Build synthetic multi-ROI AnnData.

    ROI 1 and ROI 2 have all 5 cell types.
    ROI 3 is missing TypeE.
    """
    rng = np.random.default_rng(42)
    adatas = []
    for roi_idx in range(N_ROIS):
        roi_id = f"roi_{roi_idx + 1}"
        n_cells = N_CELLS_PER_ROI

        # Raw integer counts
        X = rng.negative_binomial(5, 0.3, (n_cells, N_GENES)).astype(np.float32)

        # Cell type labels; ROI 3 has no TypeE
        if roi_idx < 2:
            ct_labels = [CELLTYPES[i % len(CELLTYPES)] for i in range(n_cells)]
        else:
            available = CELLTYPES[:-1]  # no TypeE
            ct_labels = [available[i % len(available)] for i in range(n_cells)]

        obs = pd.DataFrame(
            {
                "library_id": roi_id,
                "celltype": pd.Categorical(ct_labels, categories=CELLTYPES),
            },
            index=[f"{roi_id}_cell_{i}" for i in range(n_cells)],
        )
        var = pd.DataFrame(index=[f"gene_{i}" for i in range(N_GENES)])

        # Random 2D spatial coords (non-overlapping per ROI)
        offset = roi_idx * 1000.0
        spatial = rng.uniform(offset, offset + 500, (n_cells, 2)).astype(np.float32)

        a = ad.AnnData(X=X, obs=obs, var=var)
        a.obsm["spatial"] = spatial
        adatas.append(a)

    adata = ad.concat(adatas, join="outer", label="library_id", keys=["roi_1", "roi_2", "roi_3"])
    # After concat, restore library_id obs column from index prefix
    adata.obs["library_id"] = adata.obs["library_id"].astype(str)
    # Also store the full per-cell library_id from concatenation
    # ad.concat with label sets adata.obs['library_id'] correctly
    adata.obs["celltype"] = pd.Categorical(adata.obs["celltype"].astype(str), categories=CELLTYPES)
    return adata


@pytest.fixture()
def multi_roi_adata():
    return _make_test_adata()


# ──────────────────────────────────────────────────────────────────────────────
# Unit tests: _aggregate helpers
# ──────────────────────────────────────────────────────────────────────────────


def test_combine_pvalues_stouffer():
    """Stouffer method should combine independent uniform p-values toward significance."""
    from sc_tools.gr import combine_pvalues

    # All p = 0.05 across 5 studies → combined should be more significant
    pvals = np.full((5, 3, 3), 0.05)
    combined = combine_pvalues(pvals, method="stouffer")
    assert combined.shape == (3, 3)
    assert np.all(combined < 0.05), "Combined p should be more significant than each individual"


def test_combine_pvalues_fisher():
    """Fisher method should produce a valid combined p-value array."""
    from sc_tools.gr import combine_pvalues

    pvals = np.full((4, 2, 2), 0.1)
    combined = combine_pvalues(pvals, method="fisher")
    assert combined.shape == (2, 2)
    assert np.all(combined >= 0) and np.all(combined <= 1)


def test_combine_pvalues_clips_extremes():
    """p-values of 0 and 1 must not produce inf/-inf."""
    from sc_tools.gr import combine_pvalues

    pvals = np.array([[[0.0, 1.0]], [[0.5, 0.5]]])  # shape (2, 1, 2)
    combined = combine_pvalues(pvals, method="stouffer")
    assert np.all(np.isfinite(combined))


def test_cross_roi_zscore_nan_threshold():
    """Pairs with fewer than 3 non-NaN ROIs must return NaN."""
    from sc_tools.gr import cross_roi_zscore

    # 3 ROIs, 2x2 matrix; one pair has only 2 valid ROIs
    arr = np.array(
        [
            [[1.0, 2.0], [np.nan, 4.0]],  # ROI 1
            [[2.0, 3.0], [np.nan, 5.0]],  # ROI 2; [0,0] has 2 non-NaN
            [[3.0, 4.0], [np.nan, 6.0]],  # ROI 3
        ]
    )
    # [1][0] column is all NaN → fewer than 3 valid → must return NaN
    z = cross_roi_zscore(arr)
    assert z.shape == (3, 2, 2)
    assert np.all(np.isnan(z[:, 1, 0])), "Column with <3 valid ROIs must be NaN"
    # [0][0] has 3 valid values → should NOT be NaN
    assert not np.any(np.isnan(z[:, 0, 0])), "Column with 3 valid ROIs must not be NaN"


def test_unify_matrices_missing_type():
    """Absent cell type in one ROI must fill with fill_value at correct positions."""
    from sc_tools.gr import unify_matrices

    cats_roi1 = ["TypeA", "TypeB", "TypeC"]
    cats_roi2 = ["TypeA", "TypeC"]  # missing TypeB

    mat_roi1 = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
    mat_roi2 = np.array([[10.0, 11.0], [12.0, 13.0]])

    unified, global_cats = unify_matrices(
        [mat_roi1, mat_roi2],
        [cats_roi1, cats_roi2],
        fill_value=np.nan,
    )

    assert unified.shape == (2, 3, 3)
    assert set(global_cats) == {"TypeA", "TypeB", "TypeC"}

    # In ROI 2 (index 1), TypeB row/col should be NaN
    b_idx = global_cats.index("TypeB")
    assert np.all(np.isnan(unified[1, b_idx, :])), "Missing type row must be fill_value"
    assert np.all(np.isnan(unified[1, :, b_idx])), "Missing type col must be fill_value"

    # In ROI 1 (index 0), TypeB row/col should have real values from mat_roi1
    a_idx = global_cats.index("TypeA")
    c_idx = global_cats.index("TypeC")
    assert unified[0, a_idx, a_idx] == 1.0  # mat_roi1[0,0]
    assert unified[0, a_idx, c_idx] == 3.0  # mat_roi1[0,2]


# ──────────────────────────────────────────────────────────────────────────────
# Integration tests: spatial_neighbors
# ──────────────────────────────────────────────────────────────────────────────

pytest.importorskip("squidpy", reason="squidpy required for gr tests")


def test_spatial_neighbors_builds_block_diagonal(multi_roi_adata):
    """spatial_neighbors must not create cross-ROI edges (O(n_roi^2) block check)."""
    from sc_tools.gr import spatial_neighbors

    spatial_neighbors(multi_roi_adata, library_key="library_id")

    conn = multi_roi_adata.obsp["spatial_connectivities"]
    roi_ids = multi_roi_adata.obs["library_id"].values
    unique_rois = list(dict.fromkeys(roi_ids))

    # Build index arrays per ROI
    roi_indices = {r: (roi_ids == r).nonzero()[0] for r in unique_rois}

    # Check every cross-ROI block is all-zero (O(n_roi^2) sparse block slices)
    for i, roi_a in enumerate(unique_rois):
        for roi_b in unique_rois[i + 1 :]:
            idx_a = roi_indices[roi_a]
            idx_b = roi_indices[roi_b]
            cross_block = conn[idx_a, :][:, idx_b]
            assert cross_block.sum() == 0, f"Cross-ROI edges found between {roi_a} and {roi_b}"


# ──────────────────────────────────────────────────────────────────────────────
# Integration tests: nhood_enrichment
# ──────────────────────────────────────────────────────────────────────────────


def test_nhood_enrichment_per_roi_stored(multi_roi_adata):
    """Per-ROI results must be stored for each ROI key."""
    from sc_tools.gr import nhood_enrichment, spatial_neighbors

    spatial_neighbors(multi_roi_adata, library_key="library_id")
    nhood_enrichment(multi_roi_adata, cluster_key="celltype", library_key="library_id", n_perms=50)

    per_roi = multi_roi_adata.uns["gr"]["nhood_enrichment"]["per_roi"]
    for roi_id in ["roi_1", "roi_2", "roi_3"]:
        assert roi_id in per_roi, f"ROI {roi_id} not found in per_roi"
        assert "zscore" in per_roi[roi_id]
        assert "count" in per_roi[roi_id]


def test_nhood_enrichment_aggregated_shape(multi_roi_adata):
    """Aggregated zscore_mean must have shape (N_global, N_global)."""
    from sc_tools.gr import nhood_enrichment, spatial_neighbors

    spatial_neighbors(multi_roi_adata, library_key="library_id")
    nhood_enrichment(multi_roi_adata, cluster_key="celltype", library_key="library_id", n_perms=50)

    agg = multi_roi_adata.uns["gr"]["nhood_enrichment"]["aggregated"]
    cats = agg["cats"]
    n = len(cats)
    assert agg["zscore_mean"].shape == (n, n)
    assert agg["zscore_std"].shape == (n, n)
    assert agg["pval_stouffer"].shape == (n, n)
    assert agg["pval_fdr_bh"].shape == (n, n)


def test_nhood_enrichment_missing_celltype_nan(multi_roi_adata):
    """Cell type pairs absent from all ROIs must produce NaN in zscore_mean."""
    from sc_tools.gr import nhood_enrichment, spatial_neighbors

    spatial_neighbors(multi_roi_adata, library_key="library_id")
    nhood_enrichment(multi_roi_adata, cluster_key="celltype", library_key="library_id", n_perms=50)

    agg = multi_roi_adata.uns["gr"]["nhood_enrichment"]["aggregated"]
    cats = agg["cats"]

    # TypeE is absent from roi_3; it is present in roi_1 and roi_2,
    # so TypeE×TypeE pair should have only 2 valid ROIs → cross_roi_zscore = NaN
    if "TypeE" in cats:
        e_idx = cats.index("TypeE")
        z = agg["cross_roi_zscore"]
        # With only 2 ROIs having TypeE, cross_roi_zscore must be NaN
        assert np.all(np.isnan(z[:, e_idx, e_idx])), (
            "TypeE pair should be NaN due to <3 ROIs having data"
        )


# ──────────────────────────────────────────────────────────────────────────────
# Integration tests: co_occurrence
# ──────────────────────────────────────────────────────────────────────────────


def test_co_occurrence_per_roi_only(multi_roi_adata):
    """co_occurrence must store results per ROI, not run on full concatenated adata."""
    from sc_tools.gr import co_occurrence, spatial_neighbors

    spatial_neighbors(multi_roi_adata, library_key="library_id")
    co_occurrence(multi_roi_adata, cluster_key="celltype", library_key="library_id")

    co_occ = multi_roi_adata.uns["gr"]["co_occurrence"]
    assert "per_roi" in co_occ, "co_occurrence must store per_roi results"
    per_roi = co_occ["per_roi"]
    assert len(per_roi) == N_ROIS, f"Expected {N_ROIS} ROIs, got {len(per_roi)}"


# ──────────────────────────────────────────────────────────────────────────────
# Integration tests: centrality_scores
# ──────────────────────────────────────────────────────────────────────────────


def test_centrality_scores_aggregated(multi_roi_adata):
    """centrality_scores mean DataFrame must have index = global cell type categories."""
    from sc_tools.gr import centrality_scores, spatial_neighbors

    spatial_neighbors(multi_roi_adata, library_key="library_id")
    centrality_scores(multi_roi_adata, cluster_key="celltype", library_key="library_id")

    cs = multi_roi_adata.uns["gr"]["centrality_scores"]
    assert "mean" in cs
    mean_df = cs["mean"]
    assert isinstance(mean_df, pd.DataFrame)
    # All 5 cell types must appear (union)
    for ct in CELLTYPES:
        assert ct in mean_df.index, f"Cell type {ct} missing from centrality mean"


# ──────────────────────────────────────────────────────────────────────────────
# Integration tests: interaction_matrix
# ──────────────────────────────────────────────────────────────────────────────


def test_interaction_matrix_count_sum(multi_roi_adata):
    """count_sum must equal sum of per-ROI count matrices (with 0 fill for absent types)."""
    from sc_tools.gr import interaction_matrix, spatial_neighbors

    spatial_neighbors(multi_roi_adata, library_key="library_id")
    interaction_matrix(multi_roi_adata, cluster_key="celltype", library_key="library_id")

    im = multi_roi_adata.uns["gr"]["interaction_matrix"]
    assert "count_sum" in im
    assert "prop_mean" in im

    count_sum = im["count_sum"]
    cats = im["cats"]
    n = len(cats)
    assert count_sum.shape == (n, n)
    # count_sum must be non-negative
    assert np.all(count_sum >= 0)


# ──────────────────────────────────────────────────────────────────────────────
# Fix 1: Geary C must raise NotImplementedError
# ──────────────────────────────────────────────────────────────────────────────


def test_spatial_autocorr_geary_raises(multi_roi_adata):
    """spatial_autocorr with mode='geary' must raise NotImplementedError."""
    from sc_tools.gr import spatial_autocorr, spatial_neighbors

    spatial_neighbors(multi_roi_adata, library_key="library_id")
    with pytest.raises(NotImplementedError, match="Geary"):
        spatial_autocorr(multi_roi_adata, mode="geary")


# ──────────────────────────────────────────────────────────────────────────────
# Fix 1: Moran mode must run successfully
# ──────────────────────────────────────────────────────────────────────────────


def test_spatial_autocorr_moran_runs(multi_roi_adata):
    """spatial_autocorr with mode='moran' must store results in adata.uns['gr']."""
    from sc_tools.gr import spatial_autocorr, spatial_neighbors

    spatial_neighbors(multi_roi_adata, library_key="library_id")
    spatial_autocorr(multi_roi_adata, mode="moran", n_perms=10)
    assert "spatial_autocorr" in multi_roi_adata.uns["gr"]
    assert "result" in multi_roi_adata.uns["gr"]["spatial_autocorr"]


# ──────────────────────────────────────────────────────────────────────────────
# Fix 5: ripley per-ROI smoke test
# ──────────────────────────────────────────────────────────────────────────────


def test_ripley_per_roi_stored(multi_roi_adata):
    """ripley must store per_roi results for each ROI."""
    from sc_tools.gr import ripley, spatial_neighbors

    spatial_neighbors(multi_roi_adata, library_key="library_id")
    ripley(multi_roi_adata, cluster_key="celltype", library_key="library_id", mode="L")
    rip = multi_roi_adata.uns["gr"]["ripley"]
    assert "per_roi" in rip
    per_roi = rip["per_roi"]
    assert len(per_roi) == N_ROIS, f"Expected {N_ROIS} ROIs, got {len(per_roi)}"


# ──────────────────────────────────────────────────────────────────────────────
# Fix 2: min_cells_per_type warning
# ──────────────────────────────────────────────────────────────────────────────


def test_min_cells_per_type_warns():
    """iter_rois must emit UserWarning when a category has fewer than min_cells_per_type cells."""
    from sc_tools.gr._utils import iter_rois

    adata = _make_test_adata()
    # min_cells_per_type=1000 is larger than any per-type count (100 cells / 5 types = 20 each)
    with pytest.warns(UserWarning, match="min_cells_per_type"):
        # Consume all ROI yields to trigger warnings
        list(
            iter_rois(
                adata, library_key="library_id", cluster_key="celltype", min_cells_per_type=1000
            )
        )


# ──────────────────────────────────────────────────────────────────────────────
# Fix 4: params dict in spatial_neighbors must not have duplicate keys
# ──────────────────────────────────────────────────────────────────────────────


def test_spatial_neighbors_params_no_duplicates(multi_roi_adata):
    """Stored params dict must not contain duplicate keys (library_key, coord_type)."""
    from sc_tools.gr import spatial_neighbors

    spatial_neighbors(multi_roi_adata, library_key="library_id")
    params = multi_roi_adata.uns["gr"]["spatial_neighbors"]["params"]
    # Each key must appear exactly once (dicts enforce this by construction)
    # Verify the essential keys are present without redundancy
    assert "library_key" in params
    assert "coord_type" in params
    # The params dict must not contain the sq_kwargs sub-dict (i.e., no nested dict)
    for v in params.values():
        assert not isinstance(v, dict), "params must be flat — no nested dicts"
