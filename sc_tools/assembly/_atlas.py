"""MultiOmicAtlas class wrapping MuData with patient-centric access.

Provides:
- MultiOmicAtlas: High-level wrapper around MuData with hierarchical metadata,
  patient/sample views, and cross-modal query methods.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from anndata import AnnData

__all__ = ["MultiOmicAtlas"]


class MultiOmicAtlas:
    """Multi-omic atlas wrapping a MuData object.

    Provides patient-centric access methods, metadata management,
    and serialization to h5mu format.

    Parameters
    ----------
    mdata
        A MuData object containing multi-modal data.
    metadata
        Optional dict with 'patient_metadata' and 'sample_metadata' DataFrames.
        If None, metadata is extracted from mdata.uns.
    """

    def __init__(self, mdata, metadata: dict[str, pd.DataFrame] | None = None):
        self.mdata = mdata

        if metadata is not None:
            self.mdata.uns["patient_metadata"] = metadata.get("patient_metadata", pd.DataFrame())
            self.mdata.uns["sample_metadata"] = metadata.get("sample_metadata", pd.DataFrame())

    @classmethod
    def from_modalities(
        cls,
        modalities: dict[str, AnnData],
        *,
        subject_key: str = "subject_id",
    ) -> MultiOmicAtlas:
        """Create a MultiOmicAtlas from per-modality AnnData objects.

        Parameters
        ----------
        modalities
            Dict mapping modality name to AnnData.
        subject_key
            Column name for subject identifiers.

        Returns
        -------
        MultiOmicAtlas
        """
        from sc_tools.assembly._build import build_mudata

        mdata = build_mudata(modalities, subject_key=subject_key)
        return cls(mdata)

    @property
    def patient_metadata(self) -> pd.DataFrame:
        """Patient-level metadata DataFrame."""
        meta = self.mdata.uns.get("patient_metadata")
        if meta is None:
            return pd.DataFrame()
        if isinstance(meta, pd.DataFrame):
            return meta
        # h5mu round-trip may store as dict-of-arrays
        return pd.DataFrame(meta)

    @property
    def sample_metadata(self) -> pd.DataFrame:
        """Sample-level metadata DataFrame."""
        meta = self.mdata.uns.get("sample_metadata")
        if meta is None:
            return pd.DataFrame()
        if isinstance(meta, pd.DataFrame):
            return meta
        return pd.DataFrame(meta)

    @property
    def modalities(self) -> list[str]:
        """List of modality names."""
        return list(self.mdata.mod.keys())

    @property
    def n_obs(self) -> int:
        """Total number of cells across all modalities."""
        return sum(m.n_obs for m in self.mdata.mod.values())

    def patient_view(self, patient_id: str, *, subject_key: str = "subject_id"):
        """Subset atlas to cells from a single patient across all modalities.

        Parameters
        ----------
        patient_id
            The subject_id to filter on.
        subject_key
            Column name for subject identifiers.

        Returns
        -------
        mudata.MuData
            Subset MuData containing only cells from the specified patient.
        """
        import mudata as md

        subsets = {}
        for mod_name, mod_adata in self.mdata.mod.items():
            if subject_key in mod_adata.obs.columns:
                mask = mod_adata.obs[subject_key] == patient_id
                subsets[mod_name] = mod_adata[mask].copy()
            else:
                subsets[mod_name] = mod_adata[0:0].copy()  # empty

        result = md.MuData(subsets)
        result.update()
        return result

    def sample_view(self, sample_id: str, *, sample_key: str = "sample_id"):
        """Subset atlas to cells from a single sample across all modalities.

        Parameters
        ----------
        sample_id
            The sample_id to filter on.
        sample_key
            Column name for sample identifiers.

        Returns
        -------
        mudata.MuData
            Subset MuData containing only cells from the specified sample.
        """
        import mudata as md

        subsets = {}
        for mod_name, mod_adata in self.mdata.mod.items():
            if sample_key in mod_adata.obs.columns:
                mask = mod_adata.obs[sample_key] == sample_id
                subsets[mod_name] = mod_adata[mask].copy()
            else:
                subsets[mod_name] = mod_adata[0:0].copy()

        result = md.MuData(subsets)
        result.update()
        return result

    def embed(
        self,
        *,
        method: str = "mofa",
        n_factors: int = 15,
        **kwargs,
    ) -> np.ndarray:
        """Run joint embedding on the atlas.

        Parameters
        ----------
        method
            Embedding method name (mofa, multivi, totalvi).
        n_factors
            Number of latent factors/dimensions.
        **kwargs
            Additional arguments passed to the backend.

        Returns
        -------
        np.ndarray
            Embedding array of shape (n_cells, n_factors).
        """
        from sc_tools.assembly.embed._base import get_embedding_backend

        backend = get_embedding_backend(method)
        embedding, _meta = backend.run(self.mdata, n_factors=n_factors, **kwargs)
        return embedding

    def celltype_proportions(
        self,
        *,
        celltype_key: str = "celltype",
        group_by: str = "subject_id",
    ) -> pd.DataFrame:
        """Compute cell type proportions across modalities.

        Parameters
        ----------
        celltype_key
            Column name for cell type annotations.
        group_by
            Column name to group by (e.g. subject_id, sample_id).

        Returns
        -------
        pd.DataFrame
            Columns: [group_by, celltype_key, 'modality', 'count', 'proportion'].
        """
        from sc_tools.assembly._query import celltype_proportions

        return celltype_proportions(
            self.mdata, celltype_key=celltype_key, group_by=group_by
        )

    def save(self, path: str | Path) -> None:
        """Save atlas to h5mu format.

        Parameters
        ----------
        path
            Output file path (should end in .h5mu).
        """
        self.mdata.write(str(path))

    @classmethod
    def load(cls, path: str | Path) -> MultiOmicAtlas:
        """Load atlas from h5mu file.

        Parameters
        ----------
        path
            Path to h5mu file.

        Returns
        -------
        MultiOmicAtlas
        """
        import mudata  # lazy import per CLI-08

        mdata = mudata.read(str(path))
        return cls(mdata)
