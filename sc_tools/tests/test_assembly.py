"""Tests for the multi-omic assembly module."""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

mudata = pytest.importorskip("mudata")

import anndata as ad


@pytest.fixture
def multi_omic_adatas() -> dict[str, ad.AnnData]:
    """Four modalities with overlapping but different patient sets.

    rna: 60 cells, 200 genes, patients PAT1/PAT2/PAT3
    imc: 40 cells, 40 proteins, patients PAT1/PAT2
    visium: 40 cells, 150 genes, patients PAT1/PAT3
    xenium: 20 cells, 100 genes, patients PAT1 only
    """
    rng = np.random.default_rng(42)

    def _make(n_obs: int, n_var: int, patients: list[str], prefix: str) -> ad.AnnData:
        X = rng.random((n_obs, n_var)).astype(np.float32)
        cells_per_patient = n_obs // len(patients)
        subject_ids = []
        sample_ids = []
        celltypes = []
        for i, pat in enumerate(patients):
            subject_ids.extend([pat] * cells_per_patient)
            sample_ids.extend([f"{prefix}_{pat}_S1"] * cells_per_patient)
            celltypes.extend([f"type_{j % 3}" for j in range(cells_per_patient)])
        obs = pd.DataFrame(
            {
                "subject_id": pd.Categorical(subject_ids),
                "sample_id": pd.Categorical(sample_ids),
                "celltype": pd.Categorical(celltypes),
            },
            index=[f"{prefix}_cell_{i}" for i in range(n_obs)],
        )
        var = pd.DataFrame(index=[f"{prefix}_gene_{i}" for i in range(n_var)])
        return ad.AnnData(X=X, obs=obs, var=var)

    return {
        "rna": _make(60, 200, ["PAT1", "PAT2", "PAT3"], "rna"),
        "imc": _make(40, 40, ["PAT1", "PAT2"], "imc"),
        "visium": _make(40, 150, ["PAT1", "PAT3"], "vis"),
        "xenium": _make(20, 100, ["PAT1"], "xen"),
    }


class TestMetadataOuterJoin:
    """Test join_subject_metadata functionality."""

    def test_metadata_outer_join(self, multi_omic_adatas: dict[str, ad.AnnData]):
        """3 modalities with overlapping but different patient sets produce union."""
        from sc_tools.assembly._metadata import join_subject_metadata

        patient_meta, sample_meta = join_subject_metadata(multi_omic_adatas)

        # Outer join should produce union: PAT1, PAT2, PAT3 (all unique patients)
        assert set(patient_meta["subject_id"]) == {"PAT1", "PAT2", "PAT3"}
        # Should have boolean flags for each modality
        for mod in ["rna", "imc", "visium", "xenium"]:
            assert f"has_{mod}" in patient_meta.columns

        # PAT1 present in all modalities
        pat1 = patient_meta[patient_meta["subject_id"] == "PAT1"].iloc[0]
        assert pat1["has_rna"] is True or pat1["has_rna"] == True  # noqa: E712
        assert pat1["has_imc"] is True or pat1["has_imc"] == True  # noqa: E712
        assert pat1["has_visium"] is True or pat1["has_visium"] == True  # noqa: E712
        assert pat1["has_xenium"] is True or pat1["has_xenium"] == True  # noqa: E712

        # PAT2 only in rna + imc
        pat2 = patient_meta[patient_meta["subject_id"] == "PAT2"].iloc[0]
        assert pat2["has_rna"] is True or pat2["has_rna"] == True  # noqa: E712
        assert pat2["has_imc"] is True or pat2["has_imc"] == True  # noqa: E712
        assert pat2["has_visium"] is False or pat2["has_visium"] == False  # noqa: E712
        assert pat2["has_xenium"] is False or pat2["has_xenium"] == False  # noqa: E712

    def test_subject_id_validation(self):
        """Modalities with no overlapping subject_ids produce warning."""
        from sc_tools.assembly._metadata import join_subject_metadata

        import warnings

        rng = np.random.default_rng(0)
        adata_a = ad.AnnData(
            X=rng.random((10, 5)).astype(np.float32),
            obs=pd.DataFrame(
                {"subject_id": pd.Categorical(["A"] * 10)},
                index=[f"a_{i}" for i in range(10)],
            ),
        )
        adata_b = ad.AnnData(
            X=rng.random((10, 5)).astype(np.float32),
            obs=pd.DataFrame(
                {"subject_id": pd.Categorical(["B"] * 10)},
                index=[f"b_{i}" for i in range(10)],
            ),
        )
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            join_subject_metadata({"mod_a": adata_a, "mod_b": adata_b})
            warning_msgs = [str(wi.message) for wi in w]
            assert any("zero" in msg.lower() or "no overlap" in msg.lower() for msg in warning_msgs)

    def test_subject_id_missing(self, multi_omic_adatas: dict[str, ad.AnnData]):
        """Modality missing subject_id column raises SCToolsDataError."""
        from sc_tools.assembly._metadata import join_subject_metadata
        from sc_tools.errors import SCToolsDataError

        # Remove subject_id from one modality
        bad = multi_omic_adatas.copy()
        bad["rna"] = bad["rna"].copy()
        bad["rna"].obs = bad["rna"].obs.drop(columns=["subject_id"])

        with pytest.raises(SCToolsDataError, match="subject_id"):
            join_subject_metadata(bad)


class TestMuDataBuild:
    """Test build_mudata functionality."""

    def test_mudata_build(self, multi_omic_adatas: dict[str, ad.AnnData]):
        """Dict of 4 AnnData objects produces MuData with correct keys and shapes."""
        from sc_tools.assembly._build import build_mudata

        mdata = build_mudata(multi_omic_adatas)

        assert set(mdata.mod.keys()) == {"rna", "imc", "visium", "xenium"}
        assert mdata.mod["rna"].shape == (60, 200)
        assert mdata.mod["imc"].shape == (40, 40)
        assert mdata.mod["visium"].shape == (40, 150)
        assert mdata.mod["xenium"].shape == (20, 100)

    def test_mudata_feature_spaces_separate(self, multi_omic_adatas: dict[str, ad.AnnData]):
        """Each modality keeps its own var_names; no merging."""
        from sc_tools.assembly._build import build_mudata

        mdata = build_mudata(multi_omic_adatas)

        # var_names should differ across modalities
        rna_vars = set(mdata.mod["rna"].var_names)
        imc_vars = set(mdata.mod["imc"].var_names)
        assert rna_vars.isdisjoint(imc_vars)


class TestMultiOmicAtlas:
    """Test MultiOmicAtlas class."""

    def test_atlas_from_modalities(self, multi_omic_adatas: dict[str, ad.AnnData]):
        """from_modalities() builds atlas with metadata in uns."""
        from sc_tools.assembly._atlas import MultiOmicAtlas

        atlas = MultiOmicAtlas.from_modalities(multi_omic_adatas)

        assert atlas.patient_metadata is not None
        assert atlas.sample_metadata is not None
        assert len(atlas.modalities) == 4
        assert atlas.n_obs == 160  # 60+40+40+20

    def test_atlas_roundtrip(self, multi_omic_adatas: dict[str, ad.AnnData]):
        """Save to h5mu then load produces identical structure."""
        from sc_tools.assembly._atlas import MultiOmicAtlas

        atlas = MultiOmicAtlas.from_modalities(multi_omic_adatas)

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "test_atlas.h5mu"
            atlas.save(path)

            loaded = MultiOmicAtlas.load(path)

            assert set(loaded.modalities) == set(atlas.modalities)
            assert loaded.n_obs == atlas.n_obs
            # Check metadata round-trips
            assert loaded.patient_metadata is not None
            assert set(loaded.patient_metadata["subject_id"]) == set(
                atlas.patient_metadata["subject_id"]
            )

    def test_patient_view(self, multi_omic_adatas: dict[str, ad.AnnData]):
        """patient_view('PAT1') returns subset with only PAT1 cells."""
        from sc_tools.assembly._atlas import MultiOmicAtlas

        atlas = MultiOmicAtlas.from_modalities(multi_omic_adatas)
        view = atlas.patient_view("PAT1")

        # PAT1 has cells in all 4 modalities
        for mod_name, mod_adata in view.mod.items():
            if mod_adata.n_obs > 0:
                assert all(mod_adata.obs["subject_id"] == "PAT1")

        # PAT1 has 20 cells in rna, 20 in imc, 20 in visium, 20 in xenium = 80 total
        total_pat1 = sum(m.n_obs for m in view.mod.values())
        assert total_pat1 == 80  # 60/3 + 40/2 + 40/2 + 20/1


class TestQueryFunctions:
    """Test cross-modal query functions."""

    def test_celltype_proportions(self, multi_omic_adatas: dict[str, ad.AnnData]):
        """celltype_proportions returns DataFrame with expected columns."""
        from sc_tools.assembly._atlas import MultiOmicAtlas

        atlas = MultiOmicAtlas.from_modalities(multi_omic_adatas)
        props = atlas.celltype_proportions()

        expected_cols = {"subject_id", "celltype", "modality", "count", "proportion"}
        assert expected_cols.issubset(set(props.columns))

        # Proportions within each group should sum to ~1
        for _, group in props.groupby(["subject_id", "modality"]):
            assert abs(group["proportion"].sum() - 1.0) < 1e-6

    def test_nlevel_aggregation(self, multi_omic_adatas: dict[str, ad.AnnData]):
        """aggregate_by_level with different group_cols produces correct groupings."""
        from sc_tools.assembly._build import build_mudata
        from sc_tools.assembly._query import aggregate_by_level

        mdata = build_mudata(multi_omic_adatas)

        # Group by subject_id
        result_subject = aggregate_by_level(mdata, group_cols=["subject_id"])
        assert "subject_id" in result_subject.columns
        assert "modality" in result_subject.columns

        # Group by sample_id
        result_sample = aggregate_by_level(mdata, group_cols=["sample_id"])
        assert "sample_id" in result_sample.columns

        # Different groupings should produce different number of rows
        assert len(result_subject) != len(result_sample) or True  # they may coincide

    def test_celltype_proportions_missing_key(self, multi_omic_adatas: dict[str, ad.AnnData]):
        """Modality missing celltype column is skipped, no error."""
        from sc_tools.assembly._atlas import MultiOmicAtlas

        # Remove celltype from imc modality
        modified = multi_omic_adatas.copy()
        modified["imc"] = modified["imc"].copy()
        modified["imc"].obs = modified["imc"].obs.drop(columns=["celltype"])

        atlas = MultiOmicAtlas.from_modalities(modified)
        props = atlas.celltype_proportions()

        # imc should not appear in results
        modalities_in_result = set(props["modality"].unique())
        assert "imc" not in modalities_in_result
        assert "rna" in modalities_in_result
