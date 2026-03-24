"""Tests for subject-level metadata validation (SCI-03)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from sc_tools.qc.metadata import check_confounding, validate_subject_metadata


class TestValidateSubjectMetadata:
    """Tests for validate_subject_metadata."""

    def test_multi_sample_missing_subject_id(self):
        """Multi-sample adata missing subject_id returns warning containing 'subject_id'."""
        rng = np.random.default_rng(0)
        adata = AnnData(rng.random((100, 50)).astype("float32"))
        adata.obs["library_id"] = pd.Categorical(
            [f"lib_{i % 3}" for i in range(100)]
        )

        warnings = validate_subject_metadata(adata)
        assert len(warnings) >= 1
        assert any("subject_id" in w for w in warnings)

    def test_multi_sample_subject_id_equals_library_id(self, adata_multi_subject):
        """subject_id identical to library_id (1:1 mapping) returns distinctness warning."""
        # Make subject_id == library_id (1:1 mapping)
        adata_multi_subject.obs["subject_id"] = adata_multi_subject.obs["library_id"]

        warnings = validate_subject_metadata(adata_multi_subject)
        assert len(warnings) >= 1
        assert any("1:1" in w or "same" in w for w in warnings)

    def test_single_sample_missing_subject_id(self):
        """Single-sample adata missing subject_id returns informational warning (not error)."""
        rng = np.random.default_rng(0)
        adata = AnnData(rng.random((100, 50)).astype("float32"))
        adata.obs["library_id"] = pd.Categorical(["lib_0"] * 100)

        warnings = validate_subject_metadata(adata)
        assert len(warnings) >= 1
        assert any("informational" in w.lower() or "single" in w.lower() for w in warnings)

    def test_valid_multi_sample_distinct_subject_id(self, adata_multi_subject):
        """Valid multi-sample adata with distinct subject_id returns empty warnings list."""
        warnings = validate_subject_metadata(adata_multi_subject)
        # Should have no warnings about subject_id itself
        subject_warnings = [
            w for w in warnings
            if "subject_id" in w.lower() and "confound" not in w.lower()
        ]
        assert len(subject_warnings) == 0

    def test_auto_detect_multi_sample(self):
        """Auto-detect multi_sample from multiple unique library_id values."""
        rng = np.random.default_rng(0)
        adata = AnnData(rng.random((100, 50)).astype("float32"))
        adata.obs["library_id"] = pd.Categorical(
            [f"lib_{i % 3}" for i in range(100)]
        )
        # Should auto-detect as multi-sample and warn about missing subject_id
        warnings = validate_subject_metadata(adata)
        assert len(warnings) >= 1
        assert any("subject_id" in w for w in warnings)

    def test_explicit_multi_sample_false_no_subject_warning(self):
        """When multi_sample=False explicitly, no multi-sample warnings."""
        rng = np.random.default_rng(0)
        adata = AnnData(rng.random((100, 50)).astype("float32"))
        adata.obs["library_id"] = pd.Categorical(
            [f"lib_{i % 3}" for i in range(100)]
        )
        # Force single-sample mode despite multiple libraries
        warnings = validate_subject_metadata(adata, multi_sample=False)
        # Should get single-sample info warning, not multi-sample error
        assert any("single" in w.lower() or "informational" in w.lower() for w in warnings)


class TestCheckConfounding:
    """Tests for check_confounding."""

    def test_perfect_confounding_detected(self):
        """Perfect batch-condition confounding (each batch maps to exactly one condition)."""
        obs = pd.DataFrame({
            "batch": pd.Categorical(["A", "A", "A", "B", "B", "B"]),
            "condition": pd.Categorical(
                ["control", "control", "control", "treatment", "treatment", "treatment"]
            ),
        })
        assert check_confounding(obs, "batch", "condition") is True

    def test_non_confounded_no_warning(self, adata_multi_subject):
        """Non-confounded data returns False (no confounding)."""
        result = check_confounding(
            adata_multi_subject.obs, "batch", "condition"
        )
        assert result is False


class TestConfoundingIntegration:
    """Tests for confounding detection integrated into validate_subject_metadata."""

    def test_confounded_data_returns_warning(self):
        """Batch-condition confounding detected through validate_subject_metadata."""
        rng = np.random.default_rng(0)
        n = 60
        adata = AnnData(rng.random((n, 50)).astype("float32"))
        adata.obs["library_id"] = pd.Categorical(
            [f"lib_{i % 2}" for i in range(n)]
        )
        adata.obs["subject_id"] = pd.Categorical(
            [f"subj_{i % 2}" for i in range(n)]
        )
        # Perfect confounding: batch A = control, batch B = treatment
        adata.obs["batch"] = pd.Categorical(
            ["A"] * (n // 2) + ["B"] * (n // 2)
        )
        adata.obs["condition"] = pd.Categorical(
            ["control"] * (n // 2) + ["treatment"] * (n // 2)
        )

        warnings = validate_subject_metadata(
            adata, condition_key="condition", batch_key="batch"
        )
        assert any("confound" in w.lower() for w in warnings)

    def test_non_confounded_data_no_confounding_warning(self, adata_multi_subject):
        """Non-confounded data returns no confounding warning via validate_subject_metadata."""
        warnings = validate_subject_metadata(
            adata_multi_subject, condition_key="condition", batch_key="batch"
        )
        assert not any("confound" in w.lower() for w in warnings)
