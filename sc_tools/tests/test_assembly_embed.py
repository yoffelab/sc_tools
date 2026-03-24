"""Tests for the multi-omic assembly embedding module."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

mudata = pytest.importorskip("mudata")

import anndata as ad


@pytest.fixture
def multi_omic_adatas() -> dict[str, ad.AnnData]:
    """Four modalities with overlapping patient sets."""
    rng = np.random.default_rng(42)

    def _make(n_obs: int, n_var: int, patients: list[str], prefix: str) -> ad.AnnData:
        X = rng.random((n_obs, n_var)).astype(np.float32)
        cells_per_patient = n_obs // len(patients)
        subject_ids = []
        celltypes = []
        for i, pat in enumerate(patients):
            subject_ids.extend([pat] * cells_per_patient)
            celltypes.extend([f"type_{j % 3}" for j in range(cells_per_patient)])
        obs = pd.DataFrame(
            {
                "subject_id": pd.Categorical(subject_ids),
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


@pytest.fixture
def atlas_mdata(multi_omic_adatas):
    """Build a MuData from multi_omic_adatas."""
    from sc_tools.assembly._build import build_mudata

    return build_mudata(multi_omic_adatas)


class TestEmbeddingRegistry:
    """Test embedding backend Protocol, registry, and dispatch."""

    def test_embedding_dispatch(self):
        """register then get returns the backend."""
        from sc_tools.assembly.embed._base import (
            get_embedding_backend,
            register_embedding_backend,
        )

        class DummyBackend:
            @staticmethod
            def run(mdata, *, n_factors=15, **kwargs):
                return np.zeros((1, n_factors)), {"method": "dummy"}

        register_embedding_backend("_test_dummy", DummyBackend)
        assert get_embedding_backend("_test_dummy") is DummyBackend

    def test_invalid_method(self):
        """get_embedding_backend with unknown name raises ValueError listing available."""
        from sc_tools.assembly.embed._base import get_embedding_backend

        with pytest.raises(ValueError, match="Available"):
            get_embedding_backend("nonexistent_method_xyz")

    def test_list_methods(self):
        """list_embedding_methods returns sorted list including mofa."""
        from sc_tools.assembly.embed import list_embedding_methods

        methods = list_embedding_methods()
        assert isinstance(methods, list)
        assert "mofa" in methods
        assert methods == sorted(methods)

    def test_mofa_backend_protocol(self):
        """MofaBackend has the run static method required by Protocol."""
        from sc_tools.assembly.embed._base import EmbeddingBackend
        from sc_tools.assembly.embed._mofa import MofaBackend

        # Check structural conformance: has a callable run method
        assert hasattr(MofaBackend, "run")
        assert callable(MofaBackend.run)


class TestMofaEmbedding:
    """Test MOFA+ embedding backend integration."""

    @pytest.mark.skipif(
        not pytest.importorskip("muon", reason="muon not installed"),
        reason="muon not installed",
    )
    def test_mofa_embedding(self, atlas_mdata):
        """MofaBackend.run produces obsm['X_mofa'] with correct shape."""
        mu = pytest.importorskip("muon")
        from sc_tools.assembly.embed._mofa import MofaBackend

        n_factors = 5
        embedding, meta = MofaBackend.run(atlas_mdata, n_factors=n_factors)

        assert embedding.shape[1] == n_factors
        assert "X_mofa" in atlas_mdata.obsm
        assert meta["method"] == "mofa"
        assert meta["n_factors"] == n_factors


class TestModalityValidation:
    """Test MultiVI and TotalVI modality constraint checks."""

    def test_multivi_modality_check(self, atlas_mdata):
        """MultiviBackend.run on data without ATAC raises SCToolsDataError."""
        from sc_tools.assembly.embed._multivi import MultiviBackend
        from sc_tools.errors import SCToolsDataError

        # Our test data has rna/imc/visium/xenium -- no 'atac'
        with pytest.raises(SCToolsDataError, match="RNA and ATAC"):
            MultiviBackend.run(atlas_mdata, n_factors=10)

    def test_totalvi_modality_check(self, atlas_mdata):
        """TotalviBackend.run on data without protein raises SCToolsDataError."""
        from sc_tools.assembly.embed._totalvi import TotalviBackend
        from sc_tools.errors import SCToolsDataError

        # Our test data has no 'protein' modality
        with pytest.raises(SCToolsDataError, match="RNA and protein"):
            TotalviBackend.run(atlas_mdata, n_factors=10)


class TestAtlasEmbed:
    """Test MultiOmicAtlas.embed method."""

    @pytest.mark.skipif(
        not pytest.importorskip("muon", reason="muon not installed"),
        reason="muon not installed",
    )
    def test_atlas_embed_method(self, multi_omic_adatas):
        """atlas.embed(method='mofa') stores embedding in mdata.obsm."""
        mu = pytest.importorskip("muon")
        from sc_tools.assembly._atlas import MultiOmicAtlas

        atlas = MultiOmicAtlas.from_modalities(multi_omic_adatas)
        embedding = atlas.embed(method="mofa", n_factors=5)

        assert embedding.shape[1] == 5
        assert "X_mofa" in atlas.mdata.obsm
