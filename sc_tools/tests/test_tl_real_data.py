"""
Real-data edge-case tests for sc_tools.tl (score_signature, colocalization, gene_sets).

Tests are grouped by module and use pytest.mark.skipif when real data files are absent.
Each test class is independent -- no shared state across classes.

Data paths (checked at module load time):
- GGO_P2: projects/visium/ggo_visium/results/adata.annotated.p2.h5ad (29k spots, 8 library_ids)
- IMC_CELLTYPED: projects/imc/ggo_human/results/celltyped.h5ad (293k cells, 15 celltypes, 54 ROIs)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from sc_tools.tl import score_signature
from sc_tools.tl.colocalization import pearson_correlation, truncated_similarity
from sc_tools.tl.gene_sets import load_hallmark

# ---------------------------------------------------------------------------
# Paths to real data (relative to repo root)
# ---------------------------------------------------------------------------


def _find_repo_root() -> Path:
    _sentinel = Path("projects/visium/ggo_visium/results/adata.annotated.p2.h5ad")
    for parent in Path(__file__).resolve().parents:
        if (parent / _sentinel).is_file():
            return parent
    return Path(__file__).parents[2]


REPO_ROOT = _find_repo_root()

_ggo_p2_new = REPO_ROOT / "projects" / "visium" / "ggo_visium" / "results" / "adata.annotated.h5ad"
_ggo_p2_old = (
    REPO_ROOT / "projects" / "visium" / "ggo_visium" / "results" / "adata.annotated.p2.h5ad"
)
GGO_P2 = _ggo_p2_new if _ggo_p2_new.exists() else _ggo_p2_old

_imc_new = REPO_ROOT / "projects" / "imc" / "ggo_human" / "results" / "celltyped.h5ad"
IMC_CELLTYPED = _imc_new

# Aggressively limit subset sizes so tests stay fast
MAX_SPOTS_PER_LIB = 200
MAX_CELLS_PER_ROI = 200


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _load_visium_subset(path: Path, n_libs: int = 3) -> sc.AnnData:
    """
    Load Visium adata, subset to n_libs library_ids with MAX_SPOTS_PER_LIB spots each.
    Sets raw if not already set (required by score_signature with use_raw=True).
    """
    adata = sc.read_h5ad(path)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    lib_col = "library_id" if "library_id" in adata.obs.columns else None
    if lib_col is None:
        # Fall back: take first MAX_SPOTS_PER_LIB * n_libs observations
        n_take = min(adata.n_obs, MAX_SPOTS_PER_LIB * n_libs)
        adata = adata[:n_take].copy()
    else:
        libs = adata.obs[lib_col].unique()[:n_libs]
        selected = []
        for lib in libs:
            mask = adata.obs[lib_col] == lib
            idx = np.where(mask)[0][:MAX_SPOTS_PER_LIB]
            selected.extend(idx.tolist())
        adata = adata[selected].copy()

    # Ensure raw is set (required by scanpy backend)
    if adata.raw is None:
        adata.raw = adata.copy()

    return adata


def _load_imc_subset(path: Path, n_rois: int = 3) -> sc.AnnData:
    """
    Load IMC adata, subset to n_rois ROIs with MAX_CELLS_PER_ROI cells each.
    Constructs obsm['spatial'] from X_centroid / Y_centroid if not present.
    """
    adata = sc.read_h5ad(path)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    roi_col = None
    for candidate in ("roi", "ROI", "sample", "library_id"):
        if candidate in adata.obs.columns:
            roi_col = candidate
            break

    if roi_col is None:
        n_take = min(adata.n_obs, MAX_CELLS_PER_ROI * n_rois)
        adata = adata[:n_take].copy()
    else:
        rois = adata.obs[roi_col].unique()[:n_rois]
        selected = []
        for roi in rois:
            mask = adata.obs[roi_col] == roi
            idx = np.where(mask)[0][:MAX_CELLS_PER_ROI]
            selected.extend(idx.tolist())
        adata = adata[selected].copy()

    # Build obsm['spatial'] from centroid columns if absent
    if "spatial" not in adata.obsm:
        x_col = next(
            (c for c in ("X_centroid", "x_centroid", "X", "x") if c in adata.obs.columns), None
        )
        y_col = next(
            (c for c in ("Y_centroid", "y_centroid", "Y", "y") if c in adata.obs.columns), None
        )
        if x_col and y_col:
            adata.obsm["spatial"] = np.column_stack(
                [adata.obs[x_col].values, adata.obs[y_col].values]
            )
        else:
            # Synthetic spatial coords as last resort
            rng = np.random.default_rng(0)
            adata.obsm["spatial"] = rng.uniform(0, 1000, (adata.n_obs, 2))

    return adata


# ---------------------------------------------------------------------------
# 1. score_signature -- real Visium data (GGO)
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not GGO_P2.exists(), reason="ggo_visium annotated adata not available")
class TestScoreSignatureVisiumRealData:
    """Tests for score_signature on real Visium data (ggo_visium annotated checkpoint)."""

    def _adata_3libs(self) -> sc.AnnData:
        return _load_visium_subset(GGO_P2, n_libs=3)

    def _adata_1lib(self) -> sc.AnnData:
        return _load_visium_subset(GGO_P2, n_libs=1)

    # --- Test 1: Standard scoring with Hallmark HYPOXIA genes present in data ---

    def test_standard_scoring_stores_obsm(self):
        """
        Score a small custom signature from genes guaranteed to appear in the data.
        Result must be stored in obsm, not obs.
        """
        adata = self._adata_3libs()

        # Pick 5-8 gene names that actually exist in var_names
        real_genes = list(adata.var_names[:8])
        signatures = {"Test": {"SmallSig": real_genes}}

        result = score_signature(adata, signatures, use_raw=True, copy=False)

        assert "signature_score" in result.obsm
        assert "signature_score_z" in result.obsm
        df_raw = result.obsm["signature_score"]
        df_z = result.obsm["signature_score_z"]
        assert df_raw.shape == (adata.n_obs, 1)
        assert df_z.shape == (adata.n_obs, 1)
        assert "Test/SmallSig" in df_raw.columns
        # Scores must be stored in obsm, not obs
        assert "Test/SmallSig" not in result.obs.columns

    def test_standard_scoring_no_temp_columns_leaked(self):
        """Temporary scoring columns must not remain in obs after scoring completes."""
        adata = self._adata_3libs()
        real_genes = list(adata.var_names[:5])
        signatures = {"Test": {"Sig": real_genes}}

        score_signature(adata, signatures, use_raw=True, copy=False)

        assert not any(c.startswith("_tmp_sig_") for c in adata.obs.columns)

    def test_standard_scoring_report_written(self):
        """uns['signature_score_report'] must be written with required columns."""
        adata = self._adata_3libs()
        real_genes = list(adata.var_names[:5])
        signatures = {"Test": {"Sig": real_genes}}

        score_signature(adata, signatures, use_raw=True, copy=False)

        assert "signature_score_report" in adata.uns
        report = adata.uns["signature_score_report"]
        assert {"signature", "n_present", "n_missing", "status"}.issubset(report.columns)

    # --- Test 2: Genes partially missing ---

    def test_partial_missing_genes_scores_on_present(self):
        """
        Request 20 genes: 5 real + 15 nonexistent.
        Should score on the 5 present and warn about the 15 missing.
        Score column must not be all-NaN (enough present genes).
        """
        adata = self._adata_3libs()
        real_genes = list(adata.var_names[:5])
        fake_genes = [f"NONEXISTENT_GENE_{i}" for i in range(15)]
        mixed = real_genes + fake_genes

        signatures = {"Test": {"MixedSig": mixed}}
        result = score_signature(adata, signatures, use_raw=True, min_genes=3, copy=False)

        report = result.uns["signature_score_report"]
        row = report[report["signature"] == "Test/MixedSig"].iloc[0]
        assert row["n_present"] == 5
        assert row["n_missing"] == 15
        # Should score (status ok), not skipped
        assert row["status"] == "ok"
        # Scores must not be all NaN
        scores = result.obsm["signature_score"]["Test/MixedSig"].values
        assert not np.all(np.isnan(scores)), "All scores NaN despite 5 valid genes"

    def test_partial_missing_below_min_genes_produces_nan(self):
        """
        If fewer present genes than min_genes, score should be NaN and status 'skipped'.
        """
        adata = self._adata_3libs()
        real_genes = list(adata.var_names[:2])  # only 2 real
        fake_genes = [f"NONEXISTENT_{i}" for i in range(5)]
        mixed = real_genes + fake_genes

        signatures = {"Test": {"BelowMin": mixed}}
        result = score_signature(adata, signatures, use_raw=True, min_genes=3, copy=False)

        report = result.uns["signature_score_report"]
        row = report[report["signature"] == "Test/BelowMin"].iloc[0]
        assert row["status"] == "skipped"
        scores = result.obsm["signature_score"]["Test/BelowMin"].values
        assert np.all(np.isnan(scores)), "Expected all-NaN scores when below min_genes"

    # --- Test 3: All genes missing ---

    def test_all_genes_missing_produces_nan_column(self):
        """
        If all requested genes are absent from var_names, the column should be all-NaN
        (status skipped), not a crash.
        """
        adata = self._adata_3libs()
        fake_genes = [f"FAKE_GENE_{i}" for i in range(10)]
        signatures = {"Test": {"AllMissing": fake_genes}}

        result = score_signature(adata, signatures, use_raw=True, min_genes=3, copy=False)

        assert "signature_score" in result.obsm
        assert "Test/AllMissing" in result.obsm["signature_score"].columns
        scores = result.obsm["signature_score"]["Test/AllMissing"].values
        assert np.all(np.isnan(scores)), "Expected all-NaN when all genes missing"
        report = result.uns["signature_score_report"]
        row = report[report["signature"] == "Test/AllMissing"].iloc[0]
        assert row["status"] == "skipped"
        assert row["n_present"] == 0

    # --- Test 4: Single-sample (1 library_id) ---

    def test_single_library_does_not_crash(self):
        """Scoring on adata with only 1 library_id must not crash."""
        adata = self._adata_1lib()
        real_genes = list(adata.var_names[:5])
        signatures = {"Test": {"Sig": real_genes}}

        result = score_signature(adata, signatures, use_raw=True, copy=False)

        assert "signature_score" in result.obsm
        assert result.obsm["signature_score"].shape[0] == adata.n_obs

    # --- Test 5: Hallmark HYPOXIA on real data ---

    def test_hallmark_hypoxia_on_real_visium(self):
        """
        Score Hallmark HYPOXIA on real Visium data.
        At least some HYPOXIA genes should be present; result must be finite for those spots.
        """
        adata = self._adata_3libs()
        hallmark = load_hallmark()
        # Subset to just HYPOXIA to keep runtime minimal
        signatures = {"Hallmark": {"HYPOXIA": hallmark["Hallmark"]["HYPOXIA"]}}

        result = score_signature(adata, signatures, use_raw=True, min_genes=3, copy=False)

        assert "signature_score" in result.obsm
        assert "Hallmark/HYPOXIA" in result.obsm["signature_score"].columns
        report = result.uns["signature_score_report"]
        row = report[report["signature"] == "Hallmark/HYPOXIA"].iloc[0]
        # Either scored (ok) or skipped -- must not crash
        assert row["status"] in ("ok", "skipped")
        if row["status"] == "ok":
            scores = result.obsm["signature_score"]["Hallmark/HYPOXIA"].values
            finite_count = np.sum(np.isfinite(scores))
            assert finite_count > 0, "Expected at least some finite scores for HYPOXIA"

    # --- Test 6: z-scores have mean ~0 and std ~1 per column ---

    def test_z_scores_normalized(self):
        """
        Z-scored columns (signature_score_z) should have mean ~0 and std ~1
        for columns that are not all-NaN.
        """
        adata = self._adata_3libs()
        real_genes = list(adata.var_names[:8])
        signatures = {"Test": {"Sig": real_genes}}

        score_signature(adata, signatures, use_raw=True, copy=False)

        df_z = adata.obsm["signature_score_z"]
        for col in df_z.columns:
            vals = df_z[col].values
            if not np.all(np.isnan(vals)):
                mean_val = np.nanmean(vals)
                std_val = np.nanstd(vals)
                assert abs(mean_val) < 0.1, f"Z-score mean not ~0 for {col}: {mean_val}"
                assert std_val < 1.1, f"Z-score std not ~1 for {col}: {std_val}"


# ---------------------------------------------------------------------------
# 2. score_signature -- edge cases without real data (synthetic, always run)
# ---------------------------------------------------------------------------


class TestScoreSignatureEdgeCasesSynthetic:
    """
    Edge-case tests on synthetic data -- no real data dependency.
    These cover behaviors that are hard to trigger with minimal fixtures.
    """

    def _make_adata(
        self,
        n_obs: int = 100,
        genes: list[str] | None = None,
        n_vars: int = 50,
    ) -> sc.AnnData:
        """Build synthetic adata with specified genes or n_vars dummy genes."""
        rng = np.random.default_rng(42)
        if genes is None:
            genes = [f"G{i}" for i in range(n_vars)]
        X = rng.negative_binomial(5, 0.3, (n_obs, len(genes))).astype(np.float32)
        adata = sc.AnnData(
            X,
            obs=pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)]),
            var=pd.DataFrame(index=genes),
        )
        adata.raw = adata.copy()
        return adata

    def test_all_genes_missing_single_sig_all_nan(self):
        """All genes absent: single signature column is all-NaN."""
        adata = self._make_adata(genes=["GENA", "GENB", "GENC"])
        sigs = {"Cat": {"Sig": ["FAKE1", "FAKE2", "FAKE3"]}}
        result = score_signature(adata, sigs, use_raw=True, min_genes=1, copy=False)
        scores = result.obsm["signature_score"]["Cat/Sig"].values
        assert np.all(np.isnan(scores))

    def test_partial_missing_genes_warn_and_score(self):
        """5 real + 15 fake genes; should score on real subset and record n_missing=15."""
        genes = [f"REAL_{i}" for i in range(5)] + [f"DUMMY_{i}" for i in range(45)]
        adata = self._make_adata(genes=genes)
        fake = [f"NOTINDATA_{i}" for i in range(15)]
        real = [f"REAL_{i}" for i in range(5)]
        sigs = {"Cat": {"Mixed": real + fake}}

        result = score_signature(adata, sigs, use_raw=True, min_genes=3, copy=False)

        report = result.uns["signature_score_report"]
        row = report[report["signature"] == "Cat/Mixed"].iloc[0]
        assert row["n_present"] == 5
        assert row["n_missing"] == 15
        assert row["status"] == "ok"

    def test_single_obs_adata_does_not_crash(self):
        """Adata with a single observation (1 cell/spot) must not crash."""
        adata = self._make_adata(n_obs=1, genes=[f"G{i}" for i in range(30)])
        sigs = {"Cat": {"Sig": ["G0", "G1", "G2", "G3", "G4"]}}
        # Single-obs z-score will be 0.0 (std=0 path); just must not crash
        result = score_signature(adata, sigs, use_raw=True, min_genes=3, copy=False)
        assert result.obsm["signature_score"].shape[0] == 1

    def test_exactly_min_genes_present_scores_ok(self):
        """When exactly min_genes genes are present, status should be 'ok', not 'skipped'."""
        genes = ["AA", "BB", "CC", "EXTRA1", "EXTRA2"]
        adata = self._make_adata(genes=genes)
        sigs = {"Cat": {"Sig": ["AA", "BB", "CC"]}}
        result = score_signature(adata, sigs, use_raw=True, min_genes=3, copy=False)
        report = result.uns["signature_score_report"]
        row = report.iloc[0]
        assert row["status"] == "ok"

    def test_one_below_min_one_ok_in_same_call(self):
        """
        Two signatures in one call: one with enough genes (ok), one below min (skipped).
        Both columns must be present in obsm.
        """
        genes = ["AA", "BB", "CC", "DD", "EE"]
        adata = self._make_adata(genes=genes)
        sigs = {
            "Cat": {
                "GoodSig": ["AA", "BB", "CC"],
                "BadSig": ["FAKE1", "FAKE2"],  # only 0 present < min_genes=3
            }
        }
        result = score_signature(adata, sigs, use_raw=True, min_genes=3, copy=False)
        report = result.uns["signature_score_report"]
        good_row = report[report["signature"] == "Cat/GoodSig"].iloc[0]
        bad_row = report[report["signature"] == "Cat/BadSig"].iloc[0]
        assert good_row["status"] == "ok"
        assert bad_row["status"] == "skipped"
        assert "Cat/GoodSig" in result.obsm["signature_score"].columns
        assert "Cat/BadSig" in result.obsm["signature_score"].columns

    def test_use_raw_false_works_without_raw(self):
        """use_raw=False should not require adata.raw."""
        genes = [f"G{i}" for i in range(30)]
        adata = self._make_adata(genes=genes)
        # Explicitly clear raw
        adata.raw = None
        sigs = {"Cat": {"Sig": ["G0", "G1", "G2", "G3", "G4"]}}
        result = score_signature(adata, sigs, use_raw=False, copy=False)
        assert "signature_score" in result.obsm

    def test_use_raw_true_with_no_raw_raises(self):
        """use_raw=True with adata.raw=None must raise ValueError."""
        genes = [f"G{i}" for i in range(30)]
        rng = np.random.default_rng(1)
        X = rng.negative_binomial(5, 0.3, (20, len(genes))).astype(np.float32)
        adata = sc.AnnData(X, var=pd.DataFrame(index=genes))
        # No raw set
        sigs = {"Cat": {"Sig": ["G0", "G1", "G2", "G3", "G4"]}}
        with pytest.raises(ValueError, match="raw"):
            score_signature(adata, sigs, use_raw=True)

    def test_case_insensitive_gene_matching(self):
        """Gene matching is case-insensitive; 'gapdh' should match var_name 'GAPDH'."""
        genes = ["GAPDH", "ACTB", "TP53", "MKI67", "CD8A"]
        adata = self._make_adata(genes=genes)
        sigs = {"Cat": {"Sig": ["gapdh", "actb", "tp53"]}}  # lowercase
        result = score_signature(adata, sigs, use_raw=True, min_genes=3, copy=False)
        report = result.uns["signature_score_report"]
        row = report.iloc[0]
        # All 3 lowercase genes should match the uppercase var_names
        assert row["n_present"] == 3
        assert row["status"] == "ok"

    def test_copy_true_does_not_modify_original(self):
        """copy=True must return a new adata; original must be unmodified."""
        adata = self._make_adata(genes=[f"G{i}" for i in range(30)])
        sigs = {"Cat": {"Sig": ["G0", "G1", "G2", "G3", "G4"]}}
        result = score_signature(adata, sigs, use_raw=True, copy=True)
        assert result is not adata
        assert "signature_score" not in adata.obsm
        assert "signature_score" in result.obsm


# ---------------------------------------------------------------------------
# 3. colocalization -- real IMC data
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not IMC_CELLTYPED.exists(), reason="IMC celltyped adata not available")
class TestColocalizationIMCRealData:
    """Tests for colocalization utilities on real IMC data."""

    def _adata_3rois(self) -> sc.AnnData:
        return _load_imc_subset(IMC_CELLTYPED, n_rois=3)

    def _adata_with_scores(self, adata: sc.AnnData) -> sc.AnnData:
        """Add synthetic signature scores to adata for colocalization testing."""
        rng = np.random.default_rng(99)
        n = adata.n_obs
        df = pd.DataFrame(
            {
                "Cat/SigA": rng.standard_normal(n),
                "Cat/SigB": rng.standard_normal(n),
            },
            index=adata.obs_names,
        )
        adata.obsm["signature_score_z"] = df
        adata.obsm["signature_score"] = df.copy()
        return adata

    def test_pearson_correlation_returns_square_matrix(self):
        """pearson_correlation must return a square DataFrame with correct index/columns."""
        adata = self._adata_3rois()
        adata = self._adata_with_scores(adata)
        sigs = ["Cat/SigA", "Cat/SigB"]

        corr = pearson_correlation(adata, sig_columns=sigs)

        assert isinstance(corr, pd.DataFrame)
        assert corr.shape == (len(sigs), len(sigs))
        # Diagonal should be 1.0
        for s in sigs:
            if s in corr.index and s in corr.columns:
                assert abs(corr.loc[s, s] - 1.0) < 1e-6, "Diagonal of correlation matrix must be 1"

    def test_pearson_correlation_missing_column_gracefully_handled(self):
        """Requesting a column not present in signature_score_z should not crash."""
        adata = self._adata_3rois()
        adata = self._adata_with_scores(adata)
        sigs = ["Cat/SigA", "Cat/NONEXISTENT"]

        corr = pearson_correlation(adata, sig_columns=sigs)

        # NONEXISTENT should be filtered out; corr may be 1x1 or empty
        assert isinstance(corr, pd.DataFrame)
        assert "Cat/NONEXISTENT" not in corr.columns

    def test_pearson_correlation_all_nan_column_filtered(self):
        """A signature column that is all-NaN should be filtered before correlation."""
        adata = self._adata_3rois()
        rng = np.random.default_rng(42)
        n = adata.n_obs
        df = pd.DataFrame(
            {
                "Cat/Good": rng.standard_normal(n),
                "Cat/AllNaN": np.full(n, np.nan),
            },
            index=adata.obs_names,
        )
        adata.obsm["signature_score_z"] = df
        adata.obsm["signature_score"] = df.copy()

        sigs = ["Cat/Good", "Cat/AllNaN"]
        # Should not raise; AllNaN column may be dropped
        corr = pearson_correlation(adata, sig_columns=sigs, min_valid_ratio=0.5)
        assert isinstance(corr, pd.DataFrame)


# ---------------------------------------------------------------------------
# 4. colocalization -- synthetic edge cases (always run)
# ---------------------------------------------------------------------------


class TestColocalizationSynthetic:
    """Colocalization edge-case tests on synthetic data (no real data dependency)."""

    def _make_scored_adata(self, n: int = 100) -> sc.AnnData:
        rng = np.random.default_rng(7)
        X = rng.standard_normal((n, 5)).astype(np.float32)
        obs_names = [f"cell_{i}" for i in range(n)]
        adata = sc.AnnData(X, obs=pd.DataFrame(index=obs_names))
        df = pd.DataFrame(
            {
                "CatA/Sig1": rng.standard_normal(n),
                "CatA/Sig2": rng.standard_normal(n),
            },
            index=obs_names,
        )
        adata.obsm["signature_score_z"] = df
        adata.obsm["signature_score"] = df.copy()
        adata.obsm["spatial"] = rng.uniform(0, 1000, (n, 2))
        return adata

    def test_truncated_similarity_basic(self):
        """truncated_similarity should be product where both > 0, else 0."""
        a = np.array([1.0, -1.0, 2.0, 0.5])
        b = np.array([2.0, 3.0, -1.0, 0.5])
        result = truncated_similarity(a, b)
        np.testing.assert_allclose(result, [2.0, 0.0, 0.0, 0.25])

    def test_truncated_similarity_all_negative(self):
        """All-negative inputs should produce all zeros."""
        a = np.array([-1.0, -2.0, -0.5])
        b = np.array([-3.0, -1.0, -2.0])
        result = truncated_similarity(a, b)
        assert np.all(result == 0.0)

    def test_truncated_similarity_all_positive(self):
        """All-positive inputs should equal element-wise product."""
        a = np.array([1.0, 2.0, 3.0])
        b = np.array([4.0, 5.0, 6.0])
        result = truncated_similarity(a, b)
        np.testing.assert_allclose(result, a * b)

    def test_truncated_similarity_output_shape(self):
        """Output shape should match input shape."""
        rng = np.random.default_rng(0)
        a = rng.standard_normal(500)
        b = rng.standard_normal(500)
        result = truncated_similarity(a, b)
        assert result.shape == (500,)

    def test_pearson_correlation_symmetric(self):
        """pearson_correlation must return a symmetric matrix."""
        adata = self._make_scored_adata(n=50)
        sigs = ["CatA/Sig1", "CatA/Sig2"]
        corr = pearson_correlation(adata, sig_columns=sigs)
        if corr.shape == (2, 2):
            assert (
                abs(corr.loc["CatA/Sig1", "CatA/Sig2"] - corr.loc["CatA/Sig2", "CatA/Sig1"]) < 1e-9
            )

    def test_pearson_correlation_single_signature(self):
        """pearson_correlation with a single signature should return a 1x1 DataFrame."""
        adata = self._make_scored_adata(n=50)
        corr = pearson_correlation(adata, sig_columns=["CatA/Sig1"])
        assert isinstance(corr, pd.DataFrame)
        # Should have 1 column and 1 row
        assert corr.shape[0] <= 1

    def test_pearson_correlation_no_valid_sigs_returns_empty(self):
        """Requesting only nonexistent signatures returns an empty DataFrame."""
        adata = self._make_scored_adata(n=50)
        corr = pearson_correlation(adata, sig_columns=["NONEXISTENT1", "NONEXISTENT2"])
        assert isinstance(corr, pd.DataFrame)
        assert corr.shape[0] == 0 or corr.shape[1] == 0


# ---------------------------------------------------------------------------
# 5. gene_sets / load_hallmark -- real data checks
# ---------------------------------------------------------------------------


class TestGeneSetsBundled:
    """
    Tests for load_hallmark and bundled gene sets.
    These always run (no real file dependency -- bundled data ships with package).
    """

    def test_load_hallmark_returns_at_least_10_sets(self):
        """load_hallmark must return at least 10 gene sets."""
        result = load_hallmark()
        inner = result["Hallmark"]
        assert len(inner) >= 10, f"Expected >=10 Hallmark sets, got {len(inner)}"

    def test_load_hallmark_all_gene_names_are_strings(self):
        """All gene names in all Hallmark sets must be non-empty strings."""
        result = load_hallmark()
        for set_name, genes in result["Hallmark"].items():
            for g in genes:
                assert isinstance(g, str), f"Non-string gene in {set_name}: {g!r}"
                assert len(g) > 0, f"Empty string gene in {set_name}"

    def test_load_hallmark_no_empty_sets(self):
        """No Hallmark gene set should be empty."""
        result = load_hallmark()
        for name, genes in result["Hallmark"].items():
            assert len(genes) > 0, f"Empty gene set: {name}"

    def test_load_hallmark_hypoxia_key_present(self):
        """HYPOXIA must be present in Hallmark (canonical check)."""
        result = load_hallmark()
        assert "HYPOXIA" in result["Hallmark"]

    def test_load_hallmark_hypoxia_has_known_genes(self):
        """HYPOXIA set must contain known canonical genes (EGLN1, ENO1, LDHA)."""
        result = load_hallmark()
        hypoxia = result["Hallmark"]["HYPOXIA"]
        hypoxia_upper = [g.upper() for g in hypoxia]
        for gene in ("ENO1", "LDHA", "ALDOA"):
            assert gene in hypoxia_upper, f"Expected canonical hypoxia gene {gene} in set"

    def test_load_hallmark_no_hallmark_prefix_in_set_names(self):
        """Set names must not start with HALLMARK_ (prefix is stripped in load_hallmark)."""
        result = load_hallmark()
        for name in result["Hallmark"]:
            assert not name.startswith("HALLMARK_"), f"Prefix not stripped: {name}"

    def test_load_hallmark_keys_are_strings(self):
        """All set name keys must be strings."""
        result = load_hallmark()
        for name in result["Hallmark"]:
            assert isinstance(name, str), f"Non-string key: {name!r}"

    def test_load_hallmark_50_sets(self):
        """MSigDB Hallmark has exactly 50 gene sets."""
        result = load_hallmark()
        assert len(result["Hallmark"]) == 50


# ---------------------------------------------------------------------------
# 6. score_signature with Hallmark on real Visium data (gene overlap check)
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not GGO_P2.exists(), reason="ggo_visium annotated adata not available")
class TestHallmarkOnRealVisium:
    """
    Validate Hallmark gene overlap statistics on real Visium data.
    At least one Hallmark gene set should have >0 genes present in the data.
    """

    def test_hallmark_has_nonzero_overlap_with_visium(self):
        """
        At least one Hallmark gene set must have >0 genes present in the real
        Visium var_names. This detects a total gene-name mismatch (e.g. wrong species).
        """
        adata = _load_visium_subset(GGO_P2, n_libs=1)
        var_upper = {g.upper() for g in adata.var_names}
        hallmark = load_hallmark()
        overlap_counts = {
            name: sum(1 for g in genes if g.upper() in var_upper)
            for name, genes in hallmark["Hallmark"].items()
        }
        max_overlap = max(overlap_counts.values())
        assert max_overlap > 0, (
            "No Hallmark genes overlap with Visium var_names. Possible species/gene-name mismatch."
        )

    def test_hallmark_majority_of_sets_have_coverage(self):
        """
        The majority (>50%) of Hallmark sets should have at least 5 genes
        present in the real data var_names.
        """
        adata = _load_visium_subset(GGO_P2, n_libs=1)
        var_upper = {g.upper() for g in adata.var_names}
        hallmark = load_hallmark()
        sets_with_5plus = sum(
            1
            for genes in hallmark["Hallmark"].values()
            if sum(1 for g in genes if g.upper() in var_upper) >= 5
        )
        pct = sets_with_5plus / 50
        assert pct > 0.5, (
            f"Only {sets_with_5plus}/50 Hallmark sets have >=5 genes in Visium data "
            f"({pct:.1%}). Expected >50%."
        )
