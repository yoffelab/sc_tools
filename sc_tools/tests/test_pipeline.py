"""Unit tests for sc_tools.pipeline DAG logic."""

from __future__ import annotations

import pytest

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _reset_registry():
    """Restore the global _REGISTRY after each test."""
    from sc_tools import pipeline as pl

    original = pl._REGISTRY.copy()
    yield
    pl._REGISTRY.clear()
    pl._REGISTRY.update(original)


# ---------------------------------------------------------------------------
# PhaseSpec
# ---------------------------------------------------------------------------


class TestPhaseSpec:
    def test_defaults(self):
        from sc_tools.pipeline import PhaseSpec

        spec = PhaseSpec(label="Test Phase", depends_on=[])
        assert spec.branch == "main"
        assert spec.checkpoint is None
        assert spec.optional is False
        assert spec.iterative is False

    def test_custom_fields(self):
        from sc_tools.pipeline import PhaseSpec

        spec = PhaseSpec(
            label="Custom",
            depends_on=["qc_filter"],
            branch="custom",
            checkpoint="results/custom.h5ad",
            optional=True,
            iterative=True,
        )
        assert spec.branch == "custom"
        assert spec.checkpoint == "results/custom.h5ad"
        assert spec.optional is True
        assert spec.iterative is True


# ---------------------------------------------------------------------------
# get_dag / STANDARD_PHASES
# ---------------------------------------------------------------------------


class TestGetDag:
    def test_standard_phases_present(self):
        from sc_tools.pipeline import STANDARD_PHASES, get_dag

        dag = get_dag()
        for slug in STANDARD_PHASES:
            assert slug in dag

    def test_get_dag_returns_copy(self):
        from sc_tools.pipeline import get_dag

        dag1 = get_dag()
        dag2 = get_dag()
        assert dag1 is not dag2

    def test_all_standard_slugs(self):
        from sc_tools.pipeline import get_dag

        dag = get_dag()
        expected = {
            "ingest_raw",
            "ingest_load",
            "qc_filter",
            "metadata_attach",
            "preprocess",
            "demographics",
            "scoring",
            "celltype_manual",
            "biology",
            "meta_analysis",
        }
        assert expected.issubset(set(dag.keys()))

    def test_root_phases_have_no_deps(self):
        from sc_tools.pipeline import get_dag

        dag = get_dag()
        assert dag["ingest_raw"].depends_on == []

    def test_ingestion_branch(self):
        from sc_tools.pipeline import get_dag

        dag = get_dag()
        assert dag["ingest_raw"].branch == "ingestion"
        assert dag["ingest_load"].branch == "ingestion"

    def test_parallel_branches_from_preprocess(self):
        from sc_tools.pipeline import get_dag

        dag = get_dag()
        # demographics and scoring both depend on preprocess
        assert "preprocess" in dag["demographics"].depends_on
        assert "preprocess" in dag["scoring"].depends_on

    def test_celltype_manual_is_iterative(self):
        from sc_tools.pipeline import get_dag

        dag = get_dag()
        assert dag["celltype_manual"].iterative is True
        assert dag["celltype_manual"].optional is True

    def test_demographics_is_optional(self):
        from sc_tools.pipeline import get_dag

        dag = get_dag()
        assert dag["demographics"].optional is True

    def test_meta_analysis_is_optional(self):
        from sc_tools.pipeline import get_dag

        dag = get_dag()
        assert dag["meta_analysis"].optional is True


# ---------------------------------------------------------------------------
# get_phase
# ---------------------------------------------------------------------------


class TestGetPhase:
    def test_get_existing_phase(self):
        from sc_tools.pipeline import get_phase

        spec = get_phase("qc_filter")
        assert spec.label == "QC Filtering + Concatenation"

    def test_get_missing_phase_raises(self):
        from sc_tools.pipeline import get_phase

        with pytest.raises(KeyError, match="not registered"):
            get_phase("nonexistent_phase")


# ---------------------------------------------------------------------------
# extend_dag
# ---------------------------------------------------------------------------


class TestExtendDag:
    def test_extend_with_valid_dep(self):
        from sc_tools.pipeline import PhaseSpec, extend_dag, get_dag

        extend_dag(
            "spatial_regulon",
            PhaseSpec(
                label="Spatial Regulon Analysis",
                depends_on=["scoring"],
                branch="regulon",
                checkpoint="results/adata.regulon.h5ad",
            ),
        )
        dag = get_dag()
        assert "spatial_regulon" in dag
        assert dag["spatial_regulon"].branch == "regulon"

    def test_extend_with_invalid_dep_raises(self):
        from sc_tools.pipeline import PhaseSpec, extend_dag

        with pytest.raises(ValueError, match="Unknown depends_on"):
            extend_dag(
                "bad_phase",
                PhaseSpec(label="Bad", depends_on=["nonexistent"]),
            )

    def test_extend_root_phase(self):
        from sc_tools.pipeline import PhaseSpec, extend_dag, get_dag

        extend_dag("external_data", PhaseSpec(label="External Data", depends_on=[]))
        assert "external_data" in get_dag()

    def test_custom_phase_available_after_dep(self):
        from sc_tools.pipeline import PhaseSpec, extend_dag, get_available_next

        extend_dag(
            "regulon_analysis",
            PhaseSpec(label="Regulon", depends_on=["scoring"], branch="regulon"),
        )
        all_complete = [
            "ingest_raw",
            "ingest_load",
            "qc_filter",
            "metadata_attach",
            "preprocess",
            "scoring",
        ]
        available = get_available_next(all_complete)
        assert "regulon_analysis" in available


# ---------------------------------------------------------------------------
# get_available_next
# ---------------------------------------------------------------------------


class TestGetAvailableNext:
    def test_empty_completed_returns_roots(self):
        from sc_tools.pipeline import get_available_next

        available = get_available_next([])
        assert "ingest_raw" in available
        # No phase with unsatisfied deps should appear
        assert "qc_filter" not in available

    def test_after_ingest_raw(self):
        from sc_tools.pipeline import get_available_next

        available = get_available_next(["ingest_raw"])
        assert "ingest_load" in available
        assert "qc_filter" not in available

    def test_sequential_progression(self):
        from sc_tools.pipeline import get_available_next

        completed = ["ingest_raw", "ingest_load"]
        available = get_available_next(completed)
        assert "qc_filter" in available

    def test_branching_at_preprocess(self):
        from sc_tools.pipeline import get_available_next

        completed = [
            "ingest_raw",
            "ingest_load",
            "qc_filter",
            "metadata_attach",
            "preprocess",
        ]
        available = get_available_next(completed)
        # Both parallel branches should be available
        assert "demographics" in available
        assert "scoring" in available

    def test_non_iterative_complete_excluded(self):
        from sc_tools.pipeline import get_available_next

        # qc_filter is non-iterative; once complete it should not appear
        completed = [
            "ingest_raw",
            "ingest_load",
            "qc_filter",
            "metadata_attach",
        ]
        available = get_available_next(completed)
        assert "qc_filter" not in available

    def test_iterative_phase_re_available_after_complete(self):
        from sc_tools.pipeline import get_available_next

        completed = [
            "ingest_raw",
            "ingest_load",
            "qc_filter",
            "metadata_attach",
            "preprocess",
            "scoring",
            "celltype_manual",  # already completed
        ]
        available = get_available_next(completed)
        # celltype_manual is iterative → still available
        assert "celltype_manual" in available

    def test_all_phases_complete(self):
        from sc_tools.pipeline import get_available_next

        all_phases = [
            "ingest_raw",
            "ingest_load",
            "qc_filter",
            "metadata_attach",
            "preprocess",
            "demographics",
            "scoring",
            "celltype_manual",
            "biology",
            "meta_analysis",
        ]
        available = get_available_next(all_phases)
        # Only iterative phases should remain
        from sc_tools.pipeline import get_dag

        dag = get_dag()
        for p in available:
            assert dag[p].iterative, f"Non-iterative phase '{p}' appeared after full completion"


# ---------------------------------------------------------------------------
# get_phase_checkpoint
# ---------------------------------------------------------------------------


class TestGetPhaseCheckpoint:
    def test_simple_checkpoint(self):
        from sc_tools.pipeline import get_phase_checkpoint

        cp = get_phase_checkpoint("qc_filter")
        assert cp == "results/adata.raw.h5ad"

    def test_checkpoint_with_sample_id(self):
        from sc_tools.pipeline import get_phase_checkpoint

        cp = get_phase_checkpoint("ingest_load", sample_id="sample1")
        assert cp == "data/sample1/adata.h5ad"

    def test_none_checkpoint(self):
        from sc_tools.pipeline import get_phase_checkpoint

        cp = get_phase_checkpoint("ingest_raw")
        assert cp is None

    def test_optional_checkpoint_none(self):
        from sc_tools.pipeline import get_phase_checkpoint

        cp = get_phase_checkpoint("demographics")
        assert cp is None

    def test_scoring_checkpoint(self):
        from sc_tools.pipeline import get_phase_checkpoint

        cp = get_phase_checkpoint("scoring")
        assert cp == "results/adata.scored.h5ad"

    def test_celltype_manual_checkpoint(self):
        from sc_tools.pipeline import get_phase_checkpoint

        cp = get_phase_checkpoint("celltype_manual")
        assert cp == "results/adata.celltyped.h5ad"

    def test_missing_phase_raises(self):
        from sc_tools.pipeline import get_phase_checkpoint

        with pytest.raises(KeyError):
            get_phase_checkpoint("nonexistent")


# ---------------------------------------------------------------------------
# validate_dag
# ---------------------------------------------------------------------------


class TestValidateDag:
    def test_standard_dag_is_valid(self):
        from sc_tools.pipeline import validate_dag

        errors = validate_dag()
        assert errors == [], f"Standard DAG has errors: {errors}"

    def test_broken_dag_detected(self):
        from sc_tools.pipeline import _REGISTRY, PhaseSpec, validate_dag

        # Manually inject a bad phase (bypasses extend_dag validation)
        _REGISTRY["broken"] = PhaseSpec(label="Broken", depends_on=["nonexistent_dep"])
        errors = validate_dag()
        assert any("broken" in e for e in errors)
