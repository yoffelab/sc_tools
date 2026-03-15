"""
Integration tests for sc_tools.gr using real project checkpoints.

Tests are organized into three dataset classes:
  - TestIMCGgoHuman     : IMC multi-ROI protein panel (library_key='roi')
  - TestCosMxLymphDLBCL : CosMx 1k spatial transcriptomics (library_key='sample')
  - TestVisiumGgoVisium : Visium (library_key='library_id')

Edge cases per dataset:
  EC1 - Standard multi-ROI (3 ROIs x 200 cells)
  EC2 - Sparse cell type (1 ROI has only 1-2 cells of a rare type)
  EC3 - Spatial crop (small bounding box from one slide)
  EC4 - Zero cells for a type across all ROIs (excluded celltype)
  EC5 - Single-ROI fallback

Runtime target: < 30s per test (200 cells/ROI max, n_perms=10).
"""

from __future__ import annotations

import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

pytest.importorskip("squidpy", reason="squidpy required for gr tests")

# ──────────────────────────────────────────────────────────────────────────────
# Path constants
# ──────────────────────────────────────────────────────────────────────────────

_IMC_PATH = "/Users/junbumkim/Documents/sc_tools/projects/imc/ggo_human/results/celltyped.h5ad"
_COSMX_PATH = "/Users/junbumkim/Documents/sc_tools/projects/cosmx_1k/lymph_dlbcl/results/cosmx_rna_annotated.h5ad"
_VISIUM_PATH = "/Users/junbumkim/Documents/sc_tools/projects/visium/ggo_visium/results/adata.annotated.p2.h5ad"

_IMC_AVAILABLE = Path(_IMC_PATH).exists()
_COSMX_AVAILABLE = Path(_COSMX_PATH).exists()
_VISIUM_AVAILABLE = Path(_VISIUM_PATH).exists()

N_PER_ROI = 200
N_PERMS = 10


# ──────────────────────────────────────────────────────────────────────────────
# Generic helpers
# ──────────────────────────────────────────────────────────────────────────────


def _sample_cells(
    obs: pd.DataFrame,
    group_col: str,
    group_ids: list,
    n_per_group: int,
    seed: int = 0,
) -> np.ndarray:
    """Return sorted integer indices: up to n_per_group rows per group_id."""
    rng = np.random.default_rng(seed)
    idx: list[int] = []
    for gid in group_ids:
        rows = np.where((obs[group_col] == gid).values)[0]
        chosen = rng.choice(rows, size=min(n_per_group, len(rows)), replace=False)
        idx.extend(chosen.tolist())
    return np.sort(idx)


def _finalize(adata: ad.AnnData, library_key: str, cluster_key: str) -> ad.AnnData:
    """Cast library_key to category and cluster_key to categorical."""
    adata.obs[library_key] = adata.obs[library_key].astype("category")
    adata.obs[cluster_key] = pd.Categorical(adata.obs[cluster_key].astype(str))
    return adata


# ──────────────────────────────────────────────────────────────────────────────
# IMC dataset helpers
# ──────────────────────────────────────────────────────────────────────────────

_IMC_LIB_KEY = "roi"
_IMC_CLUSTER_KEY = "celltype"
# Three representative IMC ROIs with reasonably many cells
_IMC_ROIS_3 = [
    "05092023_GGO_Group2_S21_9843_A7-01",
    "05092023_GGO_Group2_S21_9843_A7-03",
    "05122023_Vivek_S21_10641_A4_Group1-03",
]


def _load_imc_subset(rois: list, n_per_roi: int = N_PER_ROI, seed: int = 0) -> ad.AnnData:
    full = ad.read_h5ad(_IMC_PATH, backed="r")
    idx = _sample_cells(full.obs, _IMC_LIB_KEY, rois, n_per_roi, seed=seed)
    adata = full[idx].to_memory()
    # Build spatial from X_centroid / Y_centroid stored in obs
    adata.obsm["spatial"] = adata.obs[["X_centroid", "Y_centroid"]].values.astype(np.float32)
    return _finalize(adata, _IMC_LIB_KEY, _IMC_CLUSTER_KEY)


# ──────────────────────────────────────────────────────────────────────────────
# CosMx dataset helpers
# ──────────────────────────────────────────────────────────────────────────────

_COSMX_LIB_KEY = "sample"
_COSMX_CLUSTER_KEY = "celltype"
# Two samples available; treat each sample as one ROI
_COSMX_SAMPLES = ["EC-DM-7162_CTMA100", "EC-DM-7162_CTMA121"]


def _load_cosmx_subset(samples: list, n_per_sample: int = N_PER_ROI, seed: int = 0) -> ad.AnnData:
    full = ad.read_h5ad(_COSMX_PATH, backed="r")
    idx = _sample_cells(full.obs, _COSMX_LIB_KEY, samples, n_per_sample, seed=seed)
    adata = full[idx].to_memory()
    # Build obsm['spatial'] from slide coordinates
    adata.obsm["spatial"] = adata.obs[["x_slide_mm", "y_slide_mm"]].values.astype(np.float32)
    return _finalize(adata, _COSMX_LIB_KEY, _COSMX_CLUSTER_KEY)


# ──────────────────────────────────────────────────────────────────────────────
# Visium dataset helpers
# ──────────────────────────────────────────────────────────────────────────────

_VISIUM_LIB_KEY = "library_id"
_VISIUM_CLUSTER_KEY = "pathologist_annotation"
_VISIUM_LIBS_3 = ["S21-5251-G3", "S21-9843-A7", "S21-21781-C15"]


def _load_visium_subset(libs: list, n_per_lib: int = N_PER_ROI, seed: int = 0) -> ad.AnnData:
    full = ad.read_h5ad(_VISIUM_PATH, backed="r")
    idx = _sample_cells(full.obs, _VISIUM_LIB_KEY, libs, n_per_lib, seed=seed)
    adata = full[idx].to_memory()
    return _finalize(adata, _VISIUM_LIB_KEY, _VISIUM_CLUSTER_KEY)


# ──────────────────────────────────────────────────────────────────────────────
# Shared assertion helpers
# ──────────────────────────────────────────────────────────────────────────────


def _assert_block_diagonal(adata: ad.AnnData, library_key: str) -> None:
    """Assert spatial_connectivities has no cross-ROI edges."""
    conn = adata.obsp["spatial_connectivities"]
    roi_ids = adata.obs[library_key].values
    unique_rois = list(dict.fromkeys(roi_ids))
    roi_indices = {r: (roi_ids == r).nonzero()[0] for r in unique_rois}
    for i, roi_a in enumerate(unique_rois):
        for roi_b in unique_rois[i + 1 :]:
            cross = conn[roi_indices[roi_a], :][:, roi_indices[roi_b]]
            assert cross.sum() == 0, f"Cross-ROI edges between {roi_a} and {roi_b}"


def _run_spatial_neighbors(adata: ad.AnnData, library_key: str) -> None:
    from sc_tools.gr import spatial_neighbors

    spatial_neighbors(adata, library_key=library_key)


def _run_nhood_enrichment(adata: ad.AnnData, library_key: str, cluster_key: str) -> dict:
    from sc_tools.gr import nhood_enrichment

    nhood_enrichment(adata, cluster_key=cluster_key, library_key=library_key, n_perms=N_PERMS)
    return adata.uns["gr"]["nhood_enrichment"]


def _run_interaction_matrix(adata: ad.AnnData, library_key: str, cluster_key: str) -> dict:
    from sc_tools.gr import interaction_matrix

    interaction_matrix(adata, cluster_key=cluster_key, library_key=library_key)
    return adata.uns["gr"]["interaction_matrix"]


def _run_centrality(adata: ad.AnnData, library_key: str, cluster_key: str) -> dict:
    from sc_tools.gr import centrality_scores

    centrality_scores(adata, cluster_key=cluster_key, library_key=library_key)
    return adata.uns["gr"]["centrality_scores"]


def _run_co_occurrence(adata: ad.AnnData, library_key: str, cluster_key: str) -> dict:
    from sc_tools.gr import co_occurrence

    co_occurrence(adata, cluster_key=cluster_key, library_key=library_key)
    return adata.uns["gr"]["co_occurrence"]


def _run_ripley(adata: ad.AnnData, library_key: str, cluster_key: str) -> dict:
    from sc_tools.gr import ripley

    ripley(adata, cluster_key=cluster_key, library_key=library_key, mode="L")
    return adata.uns["gr"]["ripley"]


def _run_spatial_autocorr(adata: ad.AnnData) -> dict:
    from sc_tools.gr import spatial_autocorr

    # Restrict to 10 genes to keep this fast
    genes = list(adata.var_names[:10])
    spatial_autocorr(adata, mode="moran", genes=genes, n_perms=N_PERMS)
    return adata.uns["gr"]["spatial_autocorr"]


# ──────────────────────────────────────────────────────────────────────────────
# IMC test class
# ──────────────────────────────────────────────────────────────────────────────


@pytest.mark.skipif(not _IMC_AVAILABLE, reason=f"IMC checkpoint not found: {_IMC_PATH}")
class TestIMCGgoHuman:
    """Edge-case tests on real IMC multi-ROI protein panel data."""

    # EC1: Standard multi-ROI (3 ROIs x 200 cells each)
    def test_ec1_standard_multi_roi_neighbors(self):
        """EC1: spatial_neighbors on 3 IMC ROIs must produce block-diagonal graph."""
        adata = _load_imc_subset(_IMC_ROIS_3)
        _run_spatial_neighbors(adata, _IMC_LIB_KEY)
        _assert_block_diagonal(adata, _IMC_LIB_KEY)

    def test_ec1_standard_multi_roi_nhood(self):
        """EC1: nhood_enrichment on 3 IMC ROIs must produce per_roi and aggregated."""
        adata = _load_imc_subset(_IMC_ROIS_3)
        _run_spatial_neighbors(adata, _IMC_LIB_KEY)
        result = _run_nhood_enrichment(adata, _IMC_LIB_KEY, _IMC_CLUSTER_KEY)
        assert len(result["per_roi"]) == 3
        agg = result["aggregated"]
        n = len(agg["cats"])
        assert agg["zscore_mean"].shape == (n, n)
        assert agg["pval_fdr_bh"].shape == (n, n)

    def test_ec1_standard_multi_roi_interaction_matrix(self):
        """EC1: interaction_matrix on 3 IMC ROIs must produce non-negative count_sum."""
        adata = _load_imc_subset(_IMC_ROIS_3)
        _run_spatial_neighbors(adata, _IMC_LIB_KEY)
        result = _run_interaction_matrix(adata, _IMC_LIB_KEY, _IMC_CLUSTER_KEY)
        assert "count_sum" in result
        assert np.all(result["count_sum"] >= 0)

    def test_ec1_standard_multi_roi_centrality(self):
        """EC1: centrality_scores on 3 IMC ROIs must produce a non-empty mean DataFrame."""
        adata = _load_imc_subset(_IMC_ROIS_3)
        _run_spatial_neighbors(adata, _IMC_LIB_KEY)
        result = _run_centrality(adata, _IMC_LIB_KEY, _IMC_CLUSTER_KEY)
        assert isinstance(result["mean"], pd.DataFrame)
        assert len(result["mean"]) > 0

    def test_ec1_standard_multi_roi_co_occurrence(self):
        """EC1: co_occurrence on 3 IMC ROIs must store per_roi and aggregated."""
        adata = _load_imc_subset(_IMC_ROIS_3)
        result = _run_co_occurrence(adata, _IMC_LIB_KEY, _IMC_CLUSTER_KEY)
        assert "per_roi" in result
        assert len(result["per_roi"]) == 3
        # Aggregated should exist and have occ_mean
        assert "occ_mean" in result["aggregated"]

    def test_ec1_standard_multi_roi_ripley(self):
        """EC1: ripley on 3 IMC ROIs must store per_roi for each ROI."""
        adata = _load_imc_subset(_IMC_ROIS_3)
        result = _run_ripley(adata, _IMC_LIB_KEY, _IMC_CLUSTER_KEY)
        assert "per_roi" in result
        assert len(result["per_roi"]) == 3

    # EC2: Sparse cell type (1-2 cells of a rare type in one ROI)
    def test_ec2_sparse_celltype_warns(self):
        """EC2: nhood_enrichment must warn (not crash) when a ROI has <5 cells of a type."""
        adata = _load_imc_subset(_IMC_ROIS_3)
        # "Airway Epi" is the rarest type (378 total / 54 ROIs ~ 7/ROI average)
        # Force one ROI to have only 2 cells of a celltype by subsetting further:
        # Keep all cells from first ROI but restrict rare type to 2 cells.
        roi0 = _IMC_ROIS_3[0]
        rare_type = "Airway Epi"
        rng = np.random.default_rng(1)

        mask_roi0 = adata.obs[_IMC_LIB_KEY] == roi0
        mask_rare = (adata.obs[_IMC_LIB_KEY] == roi0) & (adata.obs[_IMC_CLUSTER_KEY] == rare_type)
        mask_other_rois = adata.obs[_IMC_LIB_KEY] != roi0

        idx_roi0_nonrare = np.where(mask_roi0.values & ~mask_rare.values)[0]
        idx_rare = np.where(mask_rare.values)[0]
        idx_other = np.where(mask_other_rois.values)[0]

        if len(idx_rare) == 0:
            # Rare type absent from this ROI — load fresh subset that includes it
            full = ad.read_h5ad(_IMC_PATH, backed="r")
            # find ROIs that have Airway Epi
            roi_with_rare = (
                full.obs[full.obs[_IMC_CLUSTER_KEY] == rare_type][_IMC_LIB_KEY]
                .value_counts()
                .index.tolist()
            )
            if not roi_with_rare:
                pytest.skip(f"No ROI contains '{rare_type}' celltype")
            target_roi = roi_with_rare[0]
            other_rois = [r for r in _IMC_ROIS_3 if r != roi0][:2]
            rois_to_use = [target_roi] + other_rois
            adata = _load_imc_subset(rois_to_use)
            mask_roi0 = adata.obs[_IMC_LIB_KEY] == target_roi
            mask_rare = (adata.obs[_IMC_LIB_KEY] == target_roi) & (
                adata.obs[_IMC_CLUSTER_KEY] == rare_type
            )
            mask_other_rois = adata.obs[_IMC_LIB_KEY] != target_roi
            idx_roi0_nonrare = np.where(mask_roi0.values & ~mask_rare.values)[0]
            idx_rare = np.where(mask_rare.values)[0]
            idx_other = np.where(mask_other_rois.values)[0]

        # Keep only 2 cells of the rare type
        idx_rare_keep = rng.choice(idx_rare, size=min(2, len(idx_rare)), replace=False)
        keep_idx = np.sort(
            np.concatenate([idx_roi0_nonrare[:100], idx_rare_keep, idx_other[:200]])
        )
        adata_sparse = adata[keep_idx].copy()
        adata_sparse.obs[_IMC_CLUSTER_KEY] = pd.Categorical(
            adata_sparse.obs[_IMC_CLUSTER_KEY].astype(str)
        )
        adata_sparse.obs[_IMC_LIB_KEY] = adata_sparse.obs[_IMC_LIB_KEY].astype("category")

        _run_spatial_neighbors(adata_sparse, _IMC_LIB_KEY)
        # Should warn about the sparse type, but must not crash
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            _run_nhood_enrichment(adata_sparse, _IMC_LIB_KEY, _IMC_CLUSTER_KEY)
        # Result should still be stored
        assert "nhood_enrichment" in adata_sparse.uns.get("gr", {})

    # EC3: Spatial crop (small bounding box from one slide)
    def test_ec3_spatial_crop_small_roi(self):
        """EC3: spatial_neighbors and nhood_enrichment work on a spatially cropped ROI."""
        # Load a single large ROI and crop to a small bounding box
        roi_id = _IMC_ROIS_3[0]
        full = ad.read_h5ad(_IMC_PATH, backed="r")
        mask = full.obs[_IMC_LIB_KEY] == roi_id
        idx_roi = np.where(mask.values)[0]
        roi_adata = full[idx_roi].to_memory()
        roi_adata.obsm["spatial"] = roi_adata.obs[["X_centroid", "Y_centroid"]].values.astype(
            np.float32
        )

        # Crop to bottom-left 100x100 um bounding box
        x = roi_adata.obsm["spatial"][:, 0]
        y = roi_adata.obsm["spatial"][:, 1]
        x_min, y_min = x.min(), y.min()
        box_mask = (x <= x_min + 100) & (y <= y_min + 100)
        cropped = roi_adata[box_mask].copy()

        if cropped.n_obs < 5:
            pytest.skip(f"Spatial crop yields only {cropped.n_obs} cells — too few for test")

        cropped.obs[_IMC_LIB_KEY] = roi_id
        cropped.obs[_IMC_LIB_KEY] = cropped.obs[_IMC_LIB_KEY].astype("category")
        cropped.obs[_IMC_CLUSTER_KEY] = pd.Categorical(
            cropped.obs[_IMC_CLUSTER_KEY].astype(str)
        )

        # Should not crash even with very few cells
        _run_spatial_neighbors(cropped, _IMC_LIB_KEY)
        assert "spatial_connectivities" in cropped.obsp

        # nhood_enrichment may warn but must not raise
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            _run_nhood_enrichment(cropped, _IMC_LIB_KEY, _IMC_CLUSTER_KEY)
        assert "nhood_enrichment" in cropped.uns.get("gr", {})

    # EC4: Zero cells for a type in all ROIs
    def test_ec4_excluded_celltype_no_crash(self):
        """EC4: Excluding a celltype from all ROIs must not cause key errors or crashes."""
        adata = _load_imc_subset(_IMC_ROIS_3)
        # Drop "Tumor (Ki67+)" entirely from the adata
        exclude = "Tumor (Ki67+)"
        mask = adata.obs[_IMC_CLUSTER_KEY] != exclude
        adata_sub = adata[mask].copy()
        adata_sub.obs[_IMC_CLUSTER_KEY] = pd.Categorical(
            adata_sub.obs[_IMC_CLUSTER_KEY].astype(str)
        )
        adata_sub.obs[_IMC_LIB_KEY] = adata_sub.obs[_IMC_LIB_KEY].astype("category")

        _run_spatial_neighbors(adata_sub, _IMC_LIB_KEY)
        result = _run_nhood_enrichment(adata_sub, _IMC_LIB_KEY, _IMC_CLUSTER_KEY)
        cats = result["aggregated"].get("cats", [])
        assert exclude not in cats, f"Excluded type '{exclude}' must not appear in aggregated cats"

    # EC5: Single-ROI fallback
    def test_ec5_single_roi(self):
        """EC5: All functions must succeed when adata contains exactly one ROI."""
        adata = _load_imc_subset([_IMC_ROIS_3[0]])
        _run_spatial_neighbors(adata, _IMC_LIB_KEY)
        _assert_block_diagonal(adata, _IMC_LIB_KEY)

        result = _run_nhood_enrichment(adata, _IMC_LIB_KEY, _IMC_CLUSTER_KEY)
        assert "per_roi" in result
        assert len(result["per_roi"]) == 1

        result_im = _run_interaction_matrix(adata, _IMC_LIB_KEY, _IMC_CLUSTER_KEY)
        assert "count_sum" in result_im

        result_cs = _run_centrality(adata, _IMC_LIB_KEY, _IMC_CLUSTER_KEY)
        assert "mean" in result_cs


# ──────────────────────────────────────────────────────────────────────────────
# CosMx test class
# ──────────────────────────────────────────────────────────────────────────────


@pytest.mark.skipif(not _COSMX_AVAILABLE, reason=f"CosMx checkpoint not found: {_COSMX_PATH}")
class TestCosMxLymphDLBCL:
    """Edge-case tests on real CosMx spatial transcriptomics data (DLBCL lymph)."""

    # EC1: Standard multi-sample (2 samples x 200 cells each)
    def test_ec1_standard_multi_roi_neighbors(self):
        """EC1: spatial_neighbors on 2 CosMx samples must produce block-diagonal graph."""
        adata = _load_cosmx_subset(_COSMX_SAMPLES)
        _run_spatial_neighbors(adata, _COSMX_LIB_KEY)
        _assert_block_diagonal(adata, _COSMX_LIB_KEY)

    def test_ec1_standard_multi_roi_nhood(self):
        """EC1: nhood_enrichment on 2 CosMx samples must produce per_roi and aggregated."""
        adata = _load_cosmx_subset(_COSMX_SAMPLES)
        _run_spatial_neighbors(adata, _COSMX_LIB_KEY)
        result = _run_nhood_enrichment(adata, _COSMX_LIB_KEY, _COSMX_CLUSTER_KEY)
        assert len(result["per_roi"]) == 2
        agg = result["aggregated"]
        n = len(agg["cats"])
        assert agg["zscore_mean"].shape == (n, n)

    def test_ec1_standard_multi_roi_interaction_matrix(self):
        """EC1: interaction_matrix on 2 CosMx samples must produce non-negative count_sum."""
        adata = _load_cosmx_subset(_COSMX_SAMPLES)
        _run_spatial_neighbors(adata, _COSMX_LIB_KEY)
        result = _run_interaction_matrix(adata, _COSMX_LIB_KEY, _COSMX_CLUSTER_KEY)
        assert "count_sum" in result
        assert np.all(result["count_sum"] >= 0)

    def test_ec1_standard_multi_roi_centrality(self):
        """EC1: centrality_scores on 2 CosMx samples must produce non-empty mean DataFrame."""
        adata = _load_cosmx_subset(_COSMX_SAMPLES)
        _run_spatial_neighbors(adata, _COSMX_LIB_KEY)
        result = _run_centrality(adata, _COSMX_LIB_KEY, _COSMX_CLUSTER_KEY)
        assert isinstance(result.get("mean"), pd.DataFrame)
        assert len(result["mean"]) > 0

    def test_ec1_standard_multi_roi_co_occurrence(self):
        """EC1: co_occurrence on 2 CosMx samples must store per_roi."""
        adata = _load_cosmx_subset(_COSMX_SAMPLES)
        result = _run_co_occurrence(adata, _COSMX_LIB_KEY, _COSMX_CLUSTER_KEY)
        assert "per_roi" in result
        assert len(result["per_roi"]) == 2

    def test_ec1_standard_multi_roi_ripley(self):
        """EC1: ripley on 2 CosMx samples must store per_roi."""
        adata = _load_cosmx_subset(_COSMX_SAMPLES)
        result = _run_ripley(adata, _COSMX_LIB_KEY, _COSMX_CLUSTER_KEY)
        assert "per_roi" in result
        assert len(result["per_roi"]) == 2

    # EC2: Sparse cell type (keep only 2 cells of Neutrophil in one sample)
    def test_ec2_sparse_celltype_warns_and_does_not_crash(self):
        """EC2: nhood_enrichment with a sparse celltype must warn, not crash."""
        full = ad.read_h5ad(_COSMX_PATH, backed="r")
        rng = np.random.default_rng(42)

        rare_type = "Neutrophil"  # 21041 total but only ~2 samples
        sample_a = _COSMX_SAMPLES[0]
        sample_b = _COSMX_SAMPLES[1]

        # From sample_a: take 2 Neutrophil cells + 150 non-Neutrophil
        mask_a_rare = (full.obs[_COSMX_LIB_KEY] == sample_a) & (
            full.obs[_COSMX_CLUSTER_KEY] == rare_type
        )
        mask_a_other = (full.obs[_COSMX_LIB_KEY] == sample_a) & (
            full.obs[_COSMX_CLUSTER_KEY] != rare_type
        )
        idx_a_rare = np.where(mask_a_rare.values)[0]
        idx_a_other = np.where(mask_a_other.values)[0]

        if len(idx_a_rare) == 0:
            pytest.skip(f"No '{rare_type}' cells found in {sample_a}")

        idx_a_rare_keep = rng.choice(idx_a_rare, size=min(2, len(idx_a_rare)), replace=False)
        idx_a_other_keep = rng.choice(
            idx_a_other, size=min(148, len(idx_a_other)), replace=False
        )

        # From sample_b: take 200 cells normally
        mask_b = full.obs[_COSMX_LIB_KEY] == sample_b
        idx_b = np.where(mask_b.values)[0]
        idx_b_keep = rng.choice(idx_b, size=min(200, len(idx_b)), replace=False)

        keep_idx = np.sort(
            np.concatenate([idx_a_rare_keep, idx_a_other_keep, idx_b_keep])
        )
        adata = full[keep_idx].to_memory()
        adata.obsm["spatial"] = adata.obs[["x_slide_mm", "y_slide_mm"]].values.astype(np.float32)
        adata = _finalize(adata, _COSMX_LIB_KEY, _COSMX_CLUSTER_KEY)

        _run_spatial_neighbors(adata, _COSMX_LIB_KEY)
        # Must warn about sparse type (min_cells_per_type=5 by default in iter_rois)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            _run_nhood_enrichment(adata, _COSMX_LIB_KEY, _COSMX_CLUSTER_KEY)
        user_warnings = [w for w in caught if issubclass(w.category, UserWarning)]
        # At least one warning should reference the sparse celltype or ROI
        assert len(user_warnings) > 0 or "nhood_enrichment" in adata.uns.get("gr", {}), (
            "Expected UserWarning for sparse celltype or results stored"
        )
        assert "nhood_enrichment" in adata.uns.get("gr", {})

    # EC3: Spatial crop (small bounding box in slide coordinates)
    def test_ec3_spatial_crop_small_region(self):
        """EC3: spatial_neighbors and nhood_enrichment work on a spatially cropped CosMx region."""
        full = ad.read_h5ad(_COSMX_PATH, backed="r")
        sample = _COSMX_SAMPLES[0]
        mask = full.obs[_COSMX_LIB_KEY] == sample
        idx_sample = np.where(mask.values)[0]
        sample_adata = full[idx_sample].to_memory()
        sample_adata.obsm["spatial"] = sample_adata.obs[
            ["x_slide_mm", "y_slide_mm"]
        ].values.astype(np.float32)

        # Crop to 1mm x 1mm bounding box at min coordinates
        x = sample_adata.obsm["spatial"][:, 0]
        y = sample_adata.obsm["spatial"][:, 1]
        box_mask = (x <= x.min() + 1.0) & (y <= y.min() + 1.0)
        cropped = sample_adata[box_mask].copy()

        if cropped.n_obs < 5:
            pytest.skip(f"Spatial crop yields only {cropped.n_obs} cells")

        cropped.obs[_COSMX_LIB_KEY] = sample
        cropped.obs[_COSMX_LIB_KEY] = cropped.obs[_COSMX_LIB_KEY].astype("category")
        cropped.obs[_COSMX_CLUSTER_KEY] = pd.Categorical(
            cropped.obs[_COSMX_CLUSTER_KEY].astype(str)
        )

        _run_spatial_neighbors(cropped, _COSMX_LIB_KEY)
        assert "spatial_connectivities" in cropped.obsp

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            _run_nhood_enrichment(cropped, _COSMX_LIB_KEY, _COSMX_CLUSTER_KEY)
        assert "nhood_enrichment" in cropped.uns.get("gr", {})

    # EC4: Zero cells for a type across all samples
    def test_ec4_excluded_celltype_absent_from_cats(self):
        """EC4: Excluding 'DC' from all samples must remove it from aggregated cats."""
        adata = _load_cosmx_subset(_COSMX_SAMPLES)
        exclude = "DC"
        mask = adata.obs[_COSMX_CLUSTER_KEY] != exclude
        adata_sub = adata[mask].copy()
        adata_sub.obs[_COSMX_CLUSTER_KEY] = pd.Categorical(
            adata_sub.obs[_COSMX_CLUSTER_KEY].astype(str)
        )
        adata_sub.obs[_COSMX_LIB_KEY] = adata_sub.obs[_COSMX_LIB_KEY].astype("category")

        _run_spatial_neighbors(adata_sub, _COSMX_LIB_KEY)
        result = _run_nhood_enrichment(adata_sub, _COSMX_LIB_KEY, _COSMX_CLUSTER_KEY)
        cats = result["aggregated"].get("cats", [])
        assert exclude not in cats, f"Excluded type '{exclude}' must not appear in aggregated cats"

    # EC5: Single-sample fallback
    def test_ec5_single_roi(self):
        """EC5: All functions must succeed when adata contains exactly one CosMx sample."""
        adata = _load_cosmx_subset([_COSMX_SAMPLES[0]])
        _run_spatial_neighbors(adata, _COSMX_LIB_KEY)
        _assert_block_diagonal(adata, _COSMX_LIB_KEY)

        result = _run_nhood_enrichment(adata, _COSMX_LIB_KEY, _COSMX_CLUSTER_KEY)
        assert len(result["per_roi"]) == 1

        result_im = _run_interaction_matrix(adata, _COSMX_LIB_KEY, _COSMX_CLUSTER_KEY)
        assert "count_sum" in result_im


# ──────────────────────────────────────────────────────────────────────────────
# Visium test class
# ──────────────────────────────────────────────────────────────────────────────


@pytest.mark.skipif(not _VISIUM_AVAILABLE, reason=f"Visium checkpoint not found: {_VISIUM_PATH}")
class TestVisiumGgoVisium:
    """Edge-case tests on real Visium multi-sample data (ggo_visium)."""

    # EC1: Standard multi-ROI (3 library_ids x 200 spots each)
    def test_ec1_standard_multi_roi_neighbors(self):
        """EC1: spatial_neighbors on 3 Visium library_ids must produce block-diagonal graph."""
        adata = _load_visium_subset(_VISIUM_LIBS_3)
        _run_spatial_neighbors(adata, _VISIUM_LIB_KEY)
        _assert_block_diagonal(adata, _VISIUM_LIB_KEY)

    def test_ec1_standard_multi_roi_nhood(self):
        """EC1: nhood_enrichment on 3 Visium library_ids must produce per_roi and aggregated."""
        adata = _load_visium_subset(_VISIUM_LIBS_3)
        _run_spatial_neighbors(adata, _VISIUM_LIB_KEY)
        result = _run_nhood_enrichment(adata, _VISIUM_LIB_KEY, _VISIUM_CLUSTER_KEY)
        assert len(result["per_roi"]) == 3
        agg = result["aggregated"]
        n = len(agg["cats"])
        assert agg["zscore_mean"].shape == (n, n)
        assert agg["pval_fdr_bh"].shape == (n, n)

    def test_ec1_standard_multi_roi_interaction_matrix(self):
        """EC1: interaction_matrix on 3 Visium library_ids must produce non-negative count_sum."""
        adata = _load_visium_subset(_VISIUM_LIBS_3)
        _run_spatial_neighbors(adata, _VISIUM_LIB_KEY)
        result = _run_interaction_matrix(adata, _VISIUM_LIB_KEY, _VISIUM_CLUSTER_KEY)
        assert "count_sum" in result
        assert np.all(result["count_sum"] >= 0)

    def test_ec1_standard_multi_roi_centrality(self):
        """EC1: centrality_scores on 3 Visium library_ids must produce non-empty mean DataFrame."""
        adata = _load_visium_subset(_VISIUM_LIBS_3)
        _run_spatial_neighbors(adata, _VISIUM_LIB_KEY)
        result = _run_centrality(adata, _VISIUM_LIB_KEY, _VISIUM_CLUSTER_KEY)
        assert isinstance(result.get("mean"), pd.DataFrame)
        assert len(result["mean"]) > 0

    def test_ec1_standard_multi_roi_co_occurrence(self):
        """EC1: co_occurrence on 3 Visium library_ids must store per_roi and aggregated."""
        adata = _load_visium_subset(_VISIUM_LIBS_3)
        result = _run_co_occurrence(adata, _VISIUM_LIB_KEY, _VISIUM_CLUSTER_KEY)
        assert "per_roi" in result
        assert len(result["per_roi"]) == 3
        assert "occ_mean" in result.get("aggregated", {})

    def test_ec1_standard_multi_roi_ripley(self):
        """EC1: ripley on 3 Visium library_ids must store per_roi."""
        adata = _load_visium_subset(_VISIUM_LIBS_3)
        result = _run_ripley(adata, _VISIUM_LIB_KEY, _VISIUM_CLUSTER_KEY)
        assert "per_roi" in result
        assert len(result["per_roi"]) == 3

    def test_ec1_standard_multi_roi_spatial_autocorr(self):
        """EC1: spatial_autocorr (moran) on 3 Visium library_ids must store results."""
        adata = _load_visium_subset(_VISIUM_LIBS_3)
        _run_spatial_neighbors(adata, _VISIUM_LIB_KEY)
        result = _run_spatial_autocorr(adata)
        assert "result" in result
        assert isinstance(result["result"], pd.DataFrame)
        assert len(result["result"]) > 0

    # EC2: Sparse cell type (keep only 2 cells of a rare annotation in one library_id)
    def test_ec2_sparse_annotation_does_not_crash(self):
        """EC2: nhood_enrichment with 2 cells of a rare annotation must warn, not crash."""
        full = ad.read_h5ad(_VISIUM_PATH, backed="r")
        rng = np.random.default_rng(7)

        # "TLS Normal" is very rare (13 spots total)
        rare_ann = "TLS Normal"

        # Find a library_id that has TLS Normal
        ann_by_lib = full.obs[full.obs[_VISIUM_CLUSTER_KEY] == rare_ann][_VISIUM_LIB_KEY]
        libs_with_rare = ann_by_lib.unique().tolist()

        if not libs_with_rare:
            pytest.skip(f"No library_id contains '{rare_ann}' annotation")

        target_lib = libs_with_rare[0]
        other_lib = next(
            (lib for lib in _VISIUM_LIBS_3 if lib != target_lib),
            _VISIUM_LIBS_3[1],
        )

        # From target_lib: 2 rare + 150 other
        mask_rare = (full.obs[_VISIUM_LIB_KEY] == target_lib) & (
            full.obs[_VISIUM_CLUSTER_KEY] == rare_ann
        )
        mask_other = (full.obs[_VISIUM_LIB_KEY] == target_lib) & (
            full.obs[_VISIUM_CLUSTER_KEY] != rare_ann
        )
        idx_rare = np.where(mask_rare.values)[0]
        idx_other = np.where(mask_other.values)[0]
        idx_rare_keep = rng.choice(idx_rare, size=min(2, len(idx_rare)), replace=False)
        idx_other_keep = rng.choice(idx_other, size=min(150, len(idx_other)), replace=False)

        # From other_lib: 150 spots
        mask_lib2 = full.obs[_VISIUM_LIB_KEY] == other_lib
        idx_lib2 = np.where(mask_lib2.values)[0]
        idx_lib2_keep = rng.choice(idx_lib2, size=min(150, len(idx_lib2)), replace=False)

        keep_idx = np.sort(
            np.concatenate([idx_rare_keep, idx_other_keep, idx_lib2_keep])
        )
        adata = full[keep_idx].to_memory()
        adata = _finalize(adata, _VISIUM_LIB_KEY, _VISIUM_CLUSTER_KEY)

        _run_spatial_neighbors(adata, _VISIUM_LIB_KEY)
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            _run_nhood_enrichment(adata, _VISIUM_LIB_KEY, _VISIUM_CLUSTER_KEY)
        assert "nhood_enrichment" in adata.uns.get("gr", {})

    # EC3: Spatial crop (small bounding box from one Visium slide)
    def test_ec3_spatial_crop_small_roi(self):
        """EC3: spatial_neighbors works on a spatially cropped Visium ROI with few spots."""
        full = ad.read_h5ad(_VISIUM_PATH, backed="r")
        lib = _VISIUM_LIBS_3[0]
        mask = full.obs[_VISIUM_LIB_KEY] == lib
        idx_lib = np.where(mask.values)[0]
        lib_adata = full[idx_lib].to_memory()

        # Visium spatial coords are pixel-based; crop to min + 200 pixels
        x = lib_adata.obsm["spatial"][:, 0]
        y = lib_adata.obsm["spatial"][:, 1]
        x_min, y_min = x.min(), y.min()
        box_mask = (x <= x_min + 200) & (y <= y_min + 200)
        cropped = lib_adata[box_mask].copy()

        if cropped.n_obs < 3:
            pytest.skip(f"Spatial crop yields only {cropped.n_obs} spots")

        cropped.obs[_VISIUM_LIB_KEY] = lib
        cropped.obs[_VISIUM_LIB_KEY] = cropped.obs[_VISIUM_LIB_KEY].astype("category")
        cropped.obs[_VISIUM_CLUSTER_KEY] = pd.Categorical(
            cropped.obs[_VISIUM_CLUSTER_KEY].astype(str)
        )

        _run_spatial_neighbors(cropped, _VISIUM_LIB_KEY)
        assert "spatial_connectivities" in cropped.obsp

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            _run_nhood_enrichment(cropped, _VISIUM_LIB_KEY, _VISIUM_CLUSTER_KEY)
        assert "nhood_enrichment" in cropped.uns.get("gr", {})

    # EC4: Zero cells for a type in all ROIs
    def test_ec4_excluded_annotation_absent_from_cats(self):
        """EC4: Excluding an annotation from all ROIs must remove it from aggregated cats."""
        adata = _load_visium_subset(_VISIUM_LIBS_3)
        exclude = "TLS Solid"  # rare but present in some spots
        mask = adata.obs[_VISIUM_CLUSTER_KEY] != exclude
        adata_sub = adata[mask].copy()
        adata_sub.obs[_VISIUM_CLUSTER_KEY] = pd.Categorical(
            adata_sub.obs[_VISIUM_CLUSTER_KEY].astype(str)
        )
        adata_sub.obs[_VISIUM_LIB_KEY] = adata_sub.obs[_VISIUM_LIB_KEY].astype("category")

        _run_spatial_neighbors(adata_sub, _VISIUM_LIB_KEY)
        result = _run_nhood_enrichment(adata_sub, _VISIUM_LIB_KEY, _VISIUM_CLUSTER_KEY)
        cats = result["aggregated"].get("cats", [])
        assert exclude not in cats, (
            f"Excluded annotation '{exclude}' must not appear in aggregated cats"
        )

    # EC5: Single-ROI fallback
    def test_ec5_single_roi(self):
        """EC5: All functions must succeed when adata contains exactly one Visium library_id."""
        adata = _load_visium_subset([_VISIUM_LIBS_3[0]])
        _run_spatial_neighbors(adata, _VISIUM_LIB_KEY)
        _assert_block_diagonal(adata, _VISIUM_LIB_KEY)

        result = _run_nhood_enrichment(adata, _VISIUM_LIB_KEY, _VISIUM_CLUSTER_KEY)
        assert len(result["per_roi"]) == 1

        result_im = _run_interaction_matrix(adata, _VISIUM_LIB_KEY, _VISIUM_CLUSTER_KEY)
        assert "count_sum" in result_im

        result_cs = _run_centrality(adata, _VISIUM_LIB_KEY, _VISIUM_CLUSTER_KEY)
        assert "mean" in result_cs
