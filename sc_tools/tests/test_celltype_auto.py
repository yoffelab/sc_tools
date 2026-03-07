"""
TDD tests for sc_tools.tl.celltype automated cell typing module.

Tests are written first (Plan A TDD). All should fail at import until
the module is implemented.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------


def _make_adata(
    n_obs: int = 200, n_vars: int = 80, n_clusters: int = 5, seed: int = 42
) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    X = rng.negative_binomial(5, 0.3, (n_obs, n_vars)).astype(np.float32)
    var_names = [f"Gene{i}" for i in range(n_vars)]
    obs = pd.DataFrame(
        {
            "leiden": pd.Categorical(rng.choice([str(i) for i in range(n_clusters)], n_obs)),
            "library_id": pd.Categorical(rng.choice(["lib1", "lib2"], n_obs)),
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    return ad.AnnData(X, obs=obs, var=pd.DataFrame(index=var_names))


# ---------------------------------------------------------------------------
# TestBackendRegistry
# ---------------------------------------------------------------------------


class TestBackendRegistry:
    def test_register_and_list(self):
        from sc_tools.tl.celltype import list_celltype_methods, register_celltype_backend

        class _Dummy:
            @staticmethod
            def run(adata, *, cluster_key, store_proba, **kwargs):
                n = adata.n_obs
                labels = pd.Series(["Unknown"] * n, index=adata.obs_names)
                scores = pd.Series(np.zeros(n, dtype=np.float64), index=adata.obs_names)
                return labels, scores, None, {}

        register_celltype_backend("_test_dummy", _Dummy)
        methods = list_celltype_methods()
        assert "_test_dummy" in methods

    def test_unknown_method_raises(self):
        from sc_tools.tl.celltype._base import get_backend

        with pytest.raises(ValueError, match="Unknown"):
            get_backend("__nonexistent_xyz__")

    def test_get_backend_roundtrip(self):
        from sc_tools.tl.celltype._base import get_backend, register_celltype_backend

        class _Roundtrip:
            @staticmethod
            def run(adata, *, cluster_key, store_proba, **kwargs):
                n = adata.n_obs
                labels = pd.Series(["A"] * n, index=adata.obs_names)
                scores = pd.Series(np.ones(n, dtype=np.float64), index=adata.obs_names)
                return labels, scores, None, {}

        register_celltype_backend("_rt_test", _Roundtrip)
        backend = get_backend("_rt_test")
        assert backend is _Roundtrip


# ---------------------------------------------------------------------------
# TestStorageHelpers
# ---------------------------------------------------------------------------


class TestStorageHelpers:
    def _mock_backend_cls(self, label: str = "Macrophage", score: float = 1.0):
        """Return a backend class whose run() returns deterministic results."""

        class _Mock:
            @staticmethod
            def run(adata, *, cluster_key, store_proba, **kwargs):
                n = adata.n_obs
                labels = pd.Series([label] * n, index=adata.obs_names, dtype="object")
                scores = pd.Series(np.full(n, score, dtype=np.float64), index=adata.obs_names)
                proba_df = pd.DataFrame(
                    {label: np.full(n, score, dtype=np.float64)},
                    index=adata.obs_names,
                )
                return labels, scores, proba_df, {"test_meta": True}

        return _Mock

    def test_store_results_obs_keys(self):
        from sc_tools.tl.celltype._base import _store_results

        adata = _make_adata()
        mock_cls = self._mock_backend_cls()
        labels, scores, proba, meta = mock_cls.run(adata, cluster_key="leiden", store_proba=True)
        _store_results(adata, "mock", labels, scores, proba, meta, store_proba=False)

        assert "celltype_auto_mock" in adata.obs.columns
        assert "celltype_auto_mock_score" in adata.obs.columns
        assert adata.obs["celltype_auto_mock"].dtype.name == "category"
        assert adata.obs["celltype_auto_mock_score"].dtype == np.float64

    def test_store_results_uns_keys(self):
        from sc_tools.tl.celltype._base import _store_results

        adata = _make_adata()
        mock_cls = self._mock_backend_cls()
        labels, scores, proba, meta = mock_cls.run(adata, cluster_key="leiden", store_proba=True)
        _store_results(adata, "mock2", labels, scores, proba, meta, store_proba=False)

        assert "celltype_auto_mock2" in adata.uns
        assert adata.uns["celltype_auto_mock2"]["method"] == "mock2"
        assert "date" in adata.uns["celltype_auto_mock2"]
        assert adata.uns["celltype_auto_mock2"]["test_meta"] is True

    def test_store_proba_obsm(self):
        from sc_tools.tl.celltype._base import _store_results

        adata = _make_adata()
        mock_cls = self._mock_backend_cls()
        labels, scores, proba, meta = mock_cls.run(adata, cluster_key="leiden", store_proba=True)
        _store_results(adata, "mock3", labels, scores, proba, meta, store_proba=True)

        assert "celltype_proba_mock3" in adata.obsm

    def test_store_proba_false_no_obsm(self):
        from sc_tools.tl.celltype._base import _store_results

        adata = _make_adata()
        mock_cls = self._mock_backend_cls()
        labels, scores, proba, meta = mock_cls.run(adata, cluster_key="leiden", store_proba=True)
        _store_results(adata, "mock4", labels, scores, proba, meta, store_proba=False)

        assert "celltype_proba_mock4" not in adata.obsm


# ---------------------------------------------------------------------------
# TestDispatcher
# ---------------------------------------------------------------------------


class TestDispatcher:
    def _register_mock(self, name: str = "_dispatch_mock", label: str = "T cell"):
        from sc_tools.tl.celltype._base import register_celltype_backend

        class _M:
            @staticmethod
            def run(adata, *, cluster_key, store_proba, **kwargs):
                n = adata.n_obs
                labels = pd.Series([label] * n, index=adata.obs_names, dtype="object")
                scores = pd.Series(np.ones(n, dtype=np.float64), index=adata.obs_names)
                return labels, scores, None, {}

        register_celltype_backend(name, _M)
        return name

    def test_annotate_inplace(self):
        from sc_tools.tl.celltype import annotate_celltypes

        name = self._register_mock("_disp_inplace")
        adata = _make_adata()
        result = annotate_celltypes(adata, method=name)
        assert result is adata
        assert f"celltype_auto_{name}" in adata.obs.columns

    def test_annotate_copy(self):
        from sc_tools.tl.celltype import annotate_celltypes

        name = self._register_mock("_disp_copy")
        adata = _make_adata()
        result = annotate_celltypes(adata, method=name, copy=True)
        assert result is not adata
        assert f"celltype_auto_{name}" in result.obs.columns
        assert f"celltype_auto_{name}" not in adata.obs.columns

    def test_result_key_override(self):
        from sc_tools.tl.celltype import annotate_celltypes

        name = self._register_mock("_disp_rk")
        adata = _make_adata()
        annotate_celltypes(adata, method=name, result_key="custom_key")
        assert "celltype_auto_custom_key" in adata.obs.columns

    def test_dispatches_to_mock_backend(self):
        from sc_tools.tl.celltype import annotate_celltypes

        name = self._register_mock("_disp_val", label="NK cell")
        adata = _make_adata()
        annotate_celltypes(adata, method=name)
        assert (adata.obs[f"celltype_auto_{name}"] == "NK cell").all()


# ---------------------------------------------------------------------------
# TestScType
# ---------------------------------------------------------------------------


class TestScType:
    def _marker_db_basic(self, adata: ad.AnnData) -> dict:
        """Marker DB using first few genes in adata."""
        var = list(adata.var_names)
        return {
            "T cell": {"positive": var[:3], "negative": []},
            "B cell": {"positive": var[3:6], "negative": var[:1]},
        }

    def test_sctype_basic(self):
        from sc_tools.tl.celltype import annotate_celltypes

        adata = _make_adata()
        db = self._marker_db_basic(adata)
        annotate_celltypes(adata, method="sctype", marker_db=db)
        assert "celltype_auto_sctype" in adata.obs.columns
        assert adata.obs["celltype_auto_sctype"].notna().all()

    def test_sctype_negative_markers(self):
        from sc_tools.tl.celltype._sctype import ScTypeBackend

        adata = _make_adata()
        var = list(adata.var_names)
        # B cell has strong positive and T cell has strong positive but also used as B cell negative
        db = {
            "T cell": {"positive": var[:2], "negative": []},
            "B cell": {"positive": var[3:5], "negative": var[:2]},
        }
        labels, scores, proba, meta = ScTypeBackend.run(
            adata, cluster_key="leiden", store_proba=False, marker_db=db
        )
        assert isinstance(labels, pd.Series)
        assert len(labels) == adata.n_obs

    def test_sctype_unknown_genes(self):
        from sc_tools.tl.celltype._sctype import ScTypeBackend

        adata = _make_adata()
        # Include genes not in adata; should not crash
        db = {
            "T cell": {"positive": ["NOTEXIST1", "NOTEXIST2"], "negative": []},
            "B cell": {"positive": list(adata.var_names[:3]), "negative": []},
        }
        labels, scores, proba, meta = ScTypeBackend.run(
            adata, cluster_key="leiden", store_proba=False, marker_db=db
        )
        assert isinstance(labels, pd.Series)

    def test_sctype_assign_by_cluster(self):
        from sc_tools.tl.celltype._sctype import ScTypeBackend

        adata = _make_adata()
        db = {"T cell": {"positive": list(adata.var_names[:3]), "negative": []}}
        labels, scores, proba, meta = ScTypeBackend.run(
            adata, cluster_key="leiden", store_proba=False, marker_db=db, assign_by="cluster"
        )
        # Cells in same cluster should have same label
        for cluster in adata.obs["leiden"].unique():
            mask = adata.obs["leiden"] == cluster
            unique_labels = labels[mask].unique()
            assert len(unique_labels) == 1

    def test_sctype_assign_by_cell(self):
        from sc_tools.tl.celltype._sctype import ScTypeBackend

        adata = _make_adata(n_obs=100, n_vars=40, n_clusters=3, seed=7)
        db = {
            "T cell": {"positive": list(adata.var_names[:3]), "negative": []},
            "B cell": {"positive": list(adata.var_names[3:6]), "negative": []},
        }
        labels, scores, proba, meta = ScTypeBackend.run(
            adata, cluster_key="leiden", store_proba=False, marker_db=db, assign_by="cell"
        )
        assert len(labels) == adata.n_obs

    def test_sctype_min_score_threshold(self):
        from sc_tools.tl.celltype._sctype import ScTypeBackend

        adata = _make_adata()
        # Use a very high min_score so all cells go to Unknown
        db = {"T cell": {"positive": list(adata.var_names[:2]), "negative": []}}
        labels, scores, proba, meta = ScTypeBackend.run(
            adata,
            cluster_key="leiden",
            store_proba=False,
            marker_db=db,
            min_score=1e9,
        )
        assert (labels == "Unknown").all()

    def test_sctype_empty_markers_raises(self):
        from sc_tools.tl.celltype._sctype import ScTypeBackend

        adata = _make_adata()
        db: dict = {}
        with pytest.raises(ValueError, match="empty"):
            ScTypeBackend.run(adata, cluster_key="leiden", store_proba=False, marker_db=db)


# ---------------------------------------------------------------------------
# TestCustomGates
# ---------------------------------------------------------------------------


class TestCustomGates:
    def _protein_adata(self, n_obs: int = 100, seed: int = 0) -> ad.AnnData:
        rng = np.random.default_rng(seed)
        proteins = ["CD3", "CD4", "CD8a", "CD19", "CD56", "CD45", "PanCK", "CD31"]
        X = rng.uniform(0, 1, (n_obs, len(proteins))).astype(np.float32)
        obs = pd.DataFrame(
            {"leiden": pd.Categorical(rng.choice(["0", "1", "2"], n_obs))},
            index=[f"cell_{i}" for i in range(n_obs)],
        )
        return ad.AnnData(X, obs=obs, var=pd.DataFrame(index=proteins))

    def test_custom_gates_basic(self):
        from sc_tools.tl.celltype._custom_gates import CustomGatesBackend

        adata = self._protein_adata()
        gates = {
            "T cell": {"CD3": (0.5, None), "CD19": (None, 0.5)},
            "B cell": {"CD19": (0.5, None)},
        }
        labels, scores, proba, meta = CustomGatesBackend.run(
            adata, cluster_key="leiden", store_proba=False, gates=gates
        )
        assert isinstance(labels, pd.Series)
        assert len(labels) == adata.n_obs

    def test_custom_gates_hierarchical(self):
        from sc_tools.tl.celltype._custom_gates import CustomGatesBackend

        adata = self._protein_adata(seed=1)
        gates = {
            "T cell": {
                "CD3": (0.5, None),
                "children": {
                    "CD4 T cell": {"CD4": (0.5, None)},
                    "CD8 T cell": {"CD8a": (0.5, None)},
                },
            }
        }
        labels, scores, proba, meta = CustomGatesBackend.run(
            adata, cluster_key="leiden", store_proba=False, gates=gates
        )
        # Labels should be leaf types or parent/Unassigned
        valid = {"CD4 T cell", "CD8 T cell", "T cell", "Unassigned"}
        assert set(labels.unique()).issubset(valid | {"Unassigned"})

    def test_custom_gates_unassigned(self):
        from sc_tools.tl.celltype._custom_gates import CustomGatesBackend

        adata = self._protein_adata()
        # Impossible gate: require all CD3 > 2 (range is 0-1)
        gates = {"T cell": {"CD3": (2.0, None)}}
        labels, scores, proba, meta = CustomGatesBackend.run(
            adata, cluster_key="leiden", store_proba=False, gates=gates
        )
        assert (labels == "Unassigned").all()

    def test_custom_gates_max_expr_bound(self):
        from sc_tools.tl.celltype._custom_gates import CustomGatesBackend

        adata = self._protein_adata(seed=2)
        # Gate: CD3 must be below 0.0 (impossible) → all Unassigned
        gates = {"Non-T": {"CD3": (None, 0.0)}}
        labels, scores, proba, meta = CustomGatesBackend.run(
            adata, cluster_key="leiden", store_proba=False, gates=gates
        )
        assert (labels == "Unassigned").all()

    def test_custom_gates_normalized(self):
        from sc_tools.tl.celltype._custom_gates import CustomGatesBackend

        adata = self._protein_adata()
        gates = {"T cell": {"CD3": (0.5, None)}}
        # Should work the same whether normalize=True or False
        labels_norm, _, _, _ = CustomGatesBackend.run(
            adata, cluster_key="leiden", store_proba=False, gates=gates, normalize=True
        )
        labels_raw, _, _, _ = CustomGatesBackend.run(
            adata, cluster_key="leiden", store_proba=False, gates=gates, normalize=False
        )
        assert len(labels_norm) == len(labels_raw)


# ---------------------------------------------------------------------------
# TestCellTypistMocked
# ---------------------------------------------------------------------------


class TestCellTypistMocked:
    def test_celltypist_import_error_message(self):
        """If celltypist is not installed, ImportError must mention pip install."""
        import sys

        # Temporarily hide celltypist if present
        orig = sys.modules.get("celltypist")
        sys.modules["celltypist"] = None  # type: ignore[assignment]
        try:
            from sc_tools.tl.celltype._celltypist import CellTypistBackend

            adata = _make_adata()
            with pytest.raises(ImportError, match="pip install"):
                CellTypistBackend.run(adata, cluster_key="leiden", store_proba=False)
        finally:
            if orig is None:
                sys.modules.pop("celltypist", None)
            else:
                sys.modules["celltypist"] = orig

    def test_celltypist_skip_if_not_installed(self, monkeypatch):
        """annotate_celltypes with celltypist raises ImportError when not available."""
        import sys

        orig = sys.modules.get("celltypist")
        sys.modules["celltypist"] = None  # type: ignore[assignment]
        try:
            from sc_tools.tl.celltype import annotate_celltypes

            adata = _make_adata()
            with pytest.raises(ImportError):
                annotate_celltypes(adata, method="celltypist")
        finally:
            if orig is None:
                sys.modules.pop("celltypist", None)
            else:
                sys.modules["celltypist"] = orig


# ---------------------------------------------------------------------------
# TestApplyCelltypeMap
# ---------------------------------------------------------------------------


class TestApplyCelltypeMap:
    def _map_dict(self, n_clusters: int = 5) -> dict:
        celltypes = ["T cell", "B cell", "NK cell", "Macrophage", "Epithelial"]
        return {
            str(i): {
                "celltype": celltypes[i % len(celltypes)],
                "celltype_broad": celltypes[i % len(celltypes)].split()[0],
            }
            for i in range(n_clusters)
        }

    def test_apply_dict_map(self):
        from sc_tools.tl.celltype import apply_celltype_map

        adata = _make_adata()
        mapping = self._map_dict()
        apply_celltype_map(adata, mapping=mapping)
        assert "celltype" in adata.obs.columns
        assert "celltype_broad" in adata.obs.columns

    def test_apply_json_path(self):
        from sc_tools.tl.celltype import apply_celltype_map

        adata = _make_adata()
        mapping = self._map_dict()
        with tempfile.NamedTemporaryFile(suffix=".json", mode="w", delete=False) as f:
            json.dump(mapping, f)
            path = Path(f.name)
        try:
            apply_celltype_map(adata, mapping=path)
            assert "celltype" in adata.obs.columns
        finally:
            path.unlink(missing_ok=True)

    def test_missing_clusters_get_unknown(self):
        from sc_tools.tl.celltype import apply_celltype_map

        adata = _make_adata()
        # Only map cluster "0"
        mapping = {"0": {"celltype": "T cell", "celltype_broad": "T"}}
        apply_celltype_map(adata, mapping=mapping, unknown_label="UNKNOWN")
        # Cells not in cluster "0" should be "UNKNOWN"
        mask_not_zero = adata.obs["leiden"] != "0"
        assert (adata.obs.loc[mask_not_zero, "celltype"] == "UNKNOWN").all()

    def test_colors_set_in_uns(self):
        from sc_tools.tl.celltype import apply_celltype_map

        adata = _make_adata()
        mapping = self._map_dict()
        apply_celltype_map(adata, mapping=mapping)
        assert "celltype_colors" in adata.uns

    def test_celltype_broad_key(self):
        from sc_tools.tl.celltype import apply_celltype_map

        adata = _make_adata()
        mapping = self._map_dict()
        apply_celltype_map(adata, mapping=mapping, celltype_broad_key="lineage")
        assert "lineage" in adata.obs.columns

    def test_apply_copy_false_inplace(self):
        from sc_tools.tl.celltype import apply_celltype_map

        adata = _make_adata()
        result = apply_celltype_map(adata, mapping=self._map_dict(), copy=False)
        assert result is adata

    def test_apply_copy_true_returns_new(self):
        from sc_tools.tl.celltype import apply_celltype_map

        adata = _make_adata()
        result = apply_celltype_map(adata, mapping=self._map_dict(), copy=True)
        assert result is not adata
        assert "celltype" not in adata.obs.columns


# ---------------------------------------------------------------------------
# TestEnsemble
# ---------------------------------------------------------------------------


class TestEnsemble:
    def _setup_multi_method_adata(self) -> ad.AnnData:
        """AnnData with two pre-computed celltype_auto_* columns."""
        adata = _make_adata()
        # Simulate two prior runs
        adata.obs["celltype_auto_m1"] = pd.Categorical(
            ["T cell" if i % 2 == 0 else "B cell" for i in range(adata.n_obs)]
        )
        adata.obs["celltype_auto_m2"] = pd.Categorical(
            ["T cell" if i % 3 != 0 else "B cell" for i in range(adata.n_obs)]
        )
        return adata

    def test_ensemble_majority_vote(self):
        from sc_tools.tl.celltype import annotate_celltypes

        adata = self._setup_multi_method_adata()
        annotate_celltypes(
            adata,
            method="ensemble",
            methods=["m1", "m2"],
            strategy="majority_vote",
        )
        assert "celltype_auto_ensemble" in adata.obs.columns

    def test_ensemble_missing_method_raises(self):
        from sc_tools.tl.celltype import annotate_celltypes

        adata = _make_adata()
        with pytest.raises((KeyError, ValueError)):
            annotate_celltypes(
                adata,
                method="ensemble",
                methods=["__nonexistent_method__"],
                strategy="majority_vote",
            )


# ---------------------------------------------------------------------------
# TestListMethods
# ---------------------------------------------------------------------------


class TestListMethods:
    def test_list_returns_sorted(self):
        from sc_tools.tl.celltype import list_celltype_methods

        methods = list_celltype_methods()
        assert methods == sorted(methods)

    def test_all_registered_methods_present(self):
        from sc_tools.tl.celltype import list_celltype_methods

        methods = list_celltype_methods()
        # Core non-stub backends must be registered
        assert "sctype" in methods
        assert "celltypist" in methods
        assert "ensemble" in methods
        assert "custom_gates" in methods


# ---------------------------------------------------------------------------
# TestSchemaConsistency
# ---------------------------------------------------------------------------


class TestSchemaConsistency:
    def test_schema_keys_consistent(self):
        """After annotate_celltypes, all expected obs/uns keys present."""
        from sc_tools.tl.celltype import annotate_celltypes
        from sc_tools.tl.celltype._base import register_celltype_backend

        class _Schema:
            @staticmethod
            def run(adata, *, cluster_key, store_proba, **kwargs):
                n = adata.n_obs
                labels = pd.Series(["A"] * n, index=adata.obs_names, dtype="object")
                scores = pd.Series(np.ones(n, dtype=np.float64), index=adata.obs_names)
                proba = pd.DataFrame({"A": np.ones(n)}, index=adata.obs_names)
                return labels, scores, proba, {"extra": 42}

        register_celltype_backend("_schema_test", _Schema)
        adata = _make_adata()
        annotate_celltypes(adata, method="_schema_test", store_proba=True)

        assert "celltype_auto__schema_test" in adata.obs.columns
        assert "celltype_auto__schema_test_score" in adata.obs.columns
        assert "celltype_proba__schema_test" in adata.obsm
        assert "celltype_auto__schema_test" in adata.uns
        assert adata.uns["celltype_auto__schema_test"]["method"] == "_schema_test"
        assert adata.uns["celltype_auto__schema_test"]["extra"] == 42
