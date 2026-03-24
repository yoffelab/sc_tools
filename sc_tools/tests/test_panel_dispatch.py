"""Tests for panel-aware cell typing dispatch (SCI-04)."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from sc_tools.errors import SCToolsDataError


class TestPanelDetection:
    """Tests for panel detection and method restriction in annotate_celltypes."""

    def test_panel_celltypist_raises(self, adata_panel):
        """adata with n_vars=40 + method='celltypist' raises SCToolsDataError."""
        from sc_tools.tl.celltype.annotate import annotate_celltypes

        with pytest.raises(SCToolsDataError, match="panel detected"):
            annotate_celltypes(adata_panel, method="celltypist")

    def test_panel_sctype_proceeds(self, adata_panel):
        """adata with n_vars=40 + method='sctype' proceeds without error."""
        from sc_tools.tl.celltype.annotate import annotate_celltypes

        # Mock get_backend to avoid needing actual backend
        mock_backend = MagicMock()
        mock_backend.run.return_value = (
            pd.Series(["type_0"] * adata_panel.n_obs, index=adata_panel.obs_names),
            pd.Series([0.9] * adata_panel.n_obs, index=adata_panel.obs_names),
            None,
            {"method": "sctype"},
        )

        with patch(
            "sc_tools.tl.celltype.annotate.get_backend", return_value=mock_backend
        ):
            result = annotate_celltypes(adata_panel, method="sctype")
            assert result is adata_panel  # returns same object (copy=False)

    def test_panel_custom_gates_proceeds(self, adata_panel):
        """adata with n_vars=40 + method='custom_gates' proceeds without error."""
        from sc_tools.tl.celltype.annotate import annotate_celltypes

        mock_backend = MagicMock()
        mock_backend.run.return_value = (
            pd.Series(["type_0"] * adata_panel.n_obs, index=adata_panel.obs_names),
            pd.Series([0.9] * adata_panel.n_obs, index=adata_panel.obs_names),
            None,
            {"method": "custom_gates"},
        )

        with patch(
            "sc_tools.tl.celltype.annotate.get_backend", return_value=mock_backend
        ):
            result = annotate_celltypes(adata_panel, method="custom_gates")
            assert result is adata_panel

    def test_panel_force_method_logs_warning(self, adata_panel, caplog):
        """adata with n_vars=40 + method='celltypist' + force_method=True logs warning but does not raise."""
        from sc_tools.tl.celltype.annotate import annotate_celltypes

        mock_backend = MagicMock()
        mock_backend.run.return_value = (
            pd.Series(["type_0"] * adata_panel.n_obs, index=adata_panel.obs_names),
            pd.Series([0.9] * adata_panel.n_obs, index=adata_panel.obs_names),
            None,
            {"method": "celltypist"},
        )

        with patch(
            "sc_tools.tl.celltype.annotate.get_backend", return_value=mock_backend
        ):
            import logging

            with caplog.at_level(logging.WARNING, logger="sc_tools.tl.celltype.annotate"):
                result = annotate_celltypes(
                    adata_panel, method="celltypist", force_method=True
                )
                assert result is adata_panel
                assert any("force_method" in r.message for r in caplog.records)

    def test_whole_transcriptome_no_restriction(self):
        """adata with n_vars=18000 + method='celltypist' proceeds without error."""
        from sc_tools.tl.celltype.annotate import annotate_celltypes

        rng = np.random.default_rng(0)
        adata = AnnData(
            rng.random((50, 18000)).astype("float32"),
            obs=pd.DataFrame(
                {"leiden": pd.Categorical([str(i % 3) for i in range(50)])},
                index=[f"cell_{i}" for i in range(50)],
            ),
        )

        mock_backend = MagicMock()
        mock_backend.run.return_value = (
            pd.Series(["type_0"] * 50, index=adata.obs_names),
            pd.Series([0.9] * 50, index=adata.obs_names),
            None,
            {"method": "celltypist"},
        )

        with patch(
            "sc_tools.tl.celltype.annotate.get_backend", return_value=mock_backend
        ):
            result = annotate_celltypes(adata, method="celltypist")
            assert result is adata

    def test_panel_info_stored_in_uns(self, adata_panel):
        """panel_info dict is populated correctly in adata.uns['panel_dispatch']."""
        from sc_tools.tl.celltype.annotate import (
            PANEL_VALIDATED_METHODS,
            annotate_celltypes,
        )

        mock_backend = MagicMock()
        mock_backend.run.return_value = (
            pd.Series(["type_0"] * adata_panel.n_obs, index=adata_panel.obs_names),
            pd.Series([0.9] * adata_panel.n_obs, index=adata_panel.obs_names),
            None,
            {"method": "sctype"},
        )

        with patch(
            "sc_tools.tl.celltype.annotate.get_backend", return_value=mock_backend
        ):
            annotate_celltypes(adata_panel, method="sctype")

        assert "panel_dispatch" in adata_panel.uns
        info = adata_panel.uns["panel_dispatch"]
        assert info["panel_detected"] is True
        assert info["n_vars"] == 40
        assert info["restricted_methods"] == sorted(PANEL_VALIDATED_METHODS)

    def test_whole_transcriptome_panel_info_not_restricted(self):
        """Whole-transcriptome dataset has panel_dispatch with panel_detected=False and empty restricted_methods."""
        from sc_tools.tl.celltype.annotate import annotate_celltypes

        rng = np.random.default_rng(0)
        adata = AnnData(
            rng.random((50, 5000)).astype("float32"),
            obs=pd.DataFrame(
                {"leiden": pd.Categorical([str(i % 3) for i in range(50)])},
                index=[f"cell_{i}" for i in range(50)],
            ),
        )

        mock_backend = MagicMock()
        mock_backend.run.return_value = (
            pd.Series(["type_0"] * 50, index=adata.obs_names),
            pd.Series([0.9] * 50, index=adata.obs_names),
            None,
            {"method": "celltypist"},
        )

        with patch(
            "sc_tools.tl.celltype.annotate.get_backend", return_value=mock_backend
        ):
            annotate_celltypes(adata, method="celltypist")

        assert "panel_dispatch" in adata.uns
        info = adata.uns["panel_dispatch"]
        assert info["panel_detected"] is False
        assert info["restricted_methods"] == []

    def test_raw_n_vars_used_when_available(self):
        """When adata.raw exists with more vars, use raw.n_vars for panel detection."""
        from sc_tools.tl.celltype.annotate import annotate_celltypes

        rng = np.random.default_rng(0)
        # HVG-filtered: only 500 vars, but raw has 20000
        adata = AnnData(
            rng.random((50, 500)).astype("float32"),
            obs=pd.DataFrame(
                {"leiden": pd.Categorical([str(i % 3) for i in range(50)])},
                index=[f"cell_{i}" for i in range(50)],
            ),
        )
        # Set raw with full gene set
        adata_raw = AnnData(
            rng.random((50, 20000)).astype("float32"),
            obs=adata.obs.copy(),
        )
        adata.raw = adata_raw

        mock_backend = MagicMock()
        mock_backend.run.return_value = (
            pd.Series(["type_0"] * 50, index=adata.obs_names),
            pd.Series([0.9] * 50, index=adata.obs_names),
            None,
            {"method": "celltypist"},
        )

        # Should NOT raise even though adata.n_vars=500 < 1000,
        # because adata.raw.n_vars=20000
        with patch(
            "sc_tools.tl.celltype.annotate.get_backend", return_value=mock_backend
        ):
            result = annotate_celltypes(adata, method="celltypist")
            assert result is adata
