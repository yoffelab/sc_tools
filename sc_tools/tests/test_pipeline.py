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
        assert spec.required_obs == []
        assert spec.required_obsm == []
        assert spec.x_format == ""
        assert spec.qc_report is None
        assert spec.old_code == ""
        assert spec.human_in_loop is False

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

    def test_new_metadata_fields(self):
        from sc_tools.pipeline import PhaseSpec

        spec = PhaseSpec(
            label="Rich Phase",
            depends_on=[],
            required_obs=["sample", "batch"],
            required_obsm=["spatial", "X_pca"],
            x_format="raw counts",
            qc_report="pre_filter_qc_{date}.html",
            old_code="p1",
            human_in_loop=True,
        )
        assert spec.required_obs == ["sample", "batch"]
        assert spec.required_obsm == ["spatial", "X_pca"]
        assert spec.x_format == "raw counts"
        assert spec.qc_report == "pre_filter_qc_{date}.html"
        assert spec.old_code == "p1"
        assert spec.human_in_loop is True

    def test_standard_phases_have_old_codes(self):
        from sc_tools.pipeline import STANDARD_PHASES

        phases_with_old_code = [s for s, p in STANDARD_PHASES.items() if p.old_code]
        # All 10 standard phases should have old_code set
        assert len(phases_with_old_code) == 10

    def test_qc_filter_has_required_metadata(self):
        from sc_tools.pipeline import STANDARD_PHASES

        qc = STANDARD_PHASES["qc_filter"]
        assert "sample" in qc.required_obs
        assert "raw_data_dir" in qc.required_obs
        assert "spatial" in qc.required_obsm
        assert qc.x_format == "raw counts, concatenated"
        assert qc.qc_report is not None

    def test_human_in_loop_phases(self):
        from sc_tools.pipeline import STANDARD_PHASES

        hil_phases = [s for s, p in STANDARD_PHASES.items() if p.human_in_loop]
        assert "metadata_attach" in hil_phases
        assert "celltype_manual" in hil_phases


# ---------------------------------------------------------------------------
# get_dag / STANDARD_PHASES
# ---------------------------------------------------------------------------


class TestGetDag:
    def test_standard_phases_present(self):
        from sc_tools.pipeline import STANDARD_PHASES, get_dag

        dag = get_dag()
        for slug in STANDARD_PHASES:
            assert ("data_processing", slug) in dag

    def test_get_dag_returns_copy(self):
        from sc_tools.pipeline import get_dag

        dag1 = get_dag()
        dag2 = get_dag()
        assert dag1 is not dag2

    def test_all_standard_slugs(self):
        from sc_tools.pipeline import get_dag

        dag = get_dag()
        expected = {
            ("data_processing", "ingest_raw"),
            ("data_processing", "ingest_load"),
            ("data_processing", "qc_filter"),
            ("data_processing", "metadata_attach"),
            ("data_processing", "preprocess"),
            ("data_processing", "demographics"),
            ("data_processing", "scoring"),
            ("data_processing", "celltype_manual"),
            ("data_processing", "biology"),
            ("data_processing", "meta_analysis"),
        }
        assert expected.issubset(set(dag.keys()))

    def test_root_phases_have_no_deps(self):
        from sc_tools.pipeline import get_dag

        dag = get_dag()
        assert dag[("data_processing", "ingest_raw")].depends_on == []

    def test_ingestion_branch(self):
        from sc_tools.pipeline import get_dag

        dag = get_dag()
        assert dag[("data_processing", "ingest_raw")].branch == "ingestion"
        assert dag[("data_processing", "ingest_load")].branch == "ingestion"

    def test_parallel_branches_from_preprocess(self):
        from sc_tools.pipeline import get_dag

        dag = get_dag()
        dp = "data_processing"
        # demographics and scoring both depend on preprocess
        assert (dp, "preprocess") in dag[(dp, "demographics")].depends_on
        assert (dp, "preprocess") in dag[(dp, "scoring")].depends_on

    def test_celltype_manual_is_iterative(self):
        from sc_tools.pipeline import get_dag

        dag = get_dag()
        assert dag[("data_processing", "celltype_manual")].iterative is True
        assert dag[("data_processing", "celltype_manual")].optional is True

    def test_demographics_is_optional(self):
        from sc_tools.pipeline import get_dag

        dag = get_dag()
        assert dag[("data_processing", "demographics")].optional is True

    def test_meta_analysis_is_optional(self):
        from sc_tools.pipeline import get_dag

        dag = get_dag()
        assert dag[("data_processing", "meta_analysis")].optional is True


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
        assert ("data_processing", "spatial_regulon") in dag
        assert dag[("data_processing", "spatial_regulon")].branch == "regulon"

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
        assert ("data_processing", "external_data") in get_dag()

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
        assert ("data_processing", "regulon_analysis") in available


# ---------------------------------------------------------------------------
# get_available_next
# ---------------------------------------------------------------------------


class TestGetAvailableNext:
    """Tests for get_available_next which now returns PhaseKey tuples."""

    @staticmethod
    def _dp(slug: str) -> tuple[str, str]:
        return ("data_processing", slug)

    def test_empty_completed_returns_roots(self):
        from sc_tools.pipeline import get_available_next

        available = get_available_next([])
        assert self._dp("ingest_raw") in available
        # No phase with unsatisfied deps should appear
        assert self._dp("qc_filter") not in available

    def test_after_ingest_raw(self):
        from sc_tools.pipeline import get_available_next

        available = get_available_next(["ingest_raw"])
        assert self._dp("ingest_load") in available
        assert self._dp("qc_filter") not in available

    def test_sequential_progression(self):
        from sc_tools.pipeline import get_available_next

        completed = ["ingest_raw", "ingest_load"]
        available = get_available_next(completed)
        assert self._dp("qc_filter") in available

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
        assert self._dp("demographics") in available
        assert self._dp("scoring") in available

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
        assert self._dp("qc_filter") not in available

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
        assert self._dp("celltype_manual") in available

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
        assert cp == "results/adata.filtered.h5ad"

    def test_checkpoint_with_sample_id(self):
        from sc_tools.pipeline import get_phase_checkpoint

        cp = get_phase_checkpoint("ingest_load", sample_id="sample1")
        assert cp == "data/sample1/adata.ingested.h5ad"

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


# ---------------------------------------------------------------------------
# PhaseSpec.phase_group field
# ---------------------------------------------------------------------------


class TestPhaseSpecPhaseGroup:
    def test_default_phase_group_is_data_processing(self):
        from sc_tools.pipeline import PhaseSpec

        spec = PhaseSpec(label="Test", depends_on=[])
        assert spec.phase_group == "data_processing"

    def test_custom_phase_group(self):
        from sc_tools.pipeline import PhaseSpec

        spec = PhaseSpec(label="Discovery", depends_on=[], phase_group="discovery")
        assert spec.phase_group == "discovery"

    def test_all_standard_phases_are_data_processing(self):
        from sc_tools.pipeline import STANDARD_PHASES

        for slug, spec in STANDARD_PHASES.items():
            assert spec.phase_group == "data_processing", (
                f"Standard phase '{slug}' has phase_group='{spec.phase_group}', "
                f"expected 'data_processing'"
            )


# ---------------------------------------------------------------------------
# Tuple-keyed DAG: get_dag() returns (phase_group, subphase) keys
# ---------------------------------------------------------------------------


class TestTupleDag:
    def test_dag_keys_are_tuples(self):
        from sc_tools.pipeline import get_dag

        dag = get_dag()
        for key in dag:
            assert isinstance(key, tuple), f"DAG key {key!r} is not a tuple"
            assert len(key) == 2, f"DAG key {key!r} should have 2 elements"

    def test_all_standard_tuple_keys(self):
        from sc_tools.pipeline import get_dag

        dag = get_dag()
        expected = {
            ("data_processing", "ingest_raw"),
            ("data_processing", "ingest_load"),
            ("data_processing", "qc_filter"),
            ("data_processing", "metadata_attach"),
            ("data_processing", "preprocess"),
            ("data_processing", "demographics"),
            ("data_processing", "scoring"),
            ("data_processing", "celltype_manual"),
            ("data_processing", "biology"),
            ("data_processing", "meta_analysis"),
        }
        assert expected.issubset(set(dag.keys()))

    def test_dag_values_are_phase_specs(self):
        from sc_tools.pipeline import PhaseSpec, get_dag

        dag = get_dag()
        for key, spec in dag.items():
            assert isinstance(spec, PhaseSpec), f"Value for {key!r} is not a PhaseSpec"

    def test_depends_on_are_tuples(self):
        from sc_tools.pipeline import get_dag

        dag = get_dag()
        for key, spec in dag.items():
            for dep in spec.depends_on:
                assert isinstance(dep, tuple), f"Phase {key!r} has non-tuple dependency: {dep!r}"
                assert dep in dag, f"Phase {key!r} depends on {dep!r} which is not in the DAG"

    def test_root_phase_no_deps(self):
        from sc_tools.pipeline import get_dag

        dag = get_dag()
        assert dag[("data_processing", "ingest_raw")].depends_on == []

    def test_dependency_chain(self):
        from sc_tools.pipeline import get_dag

        dag = get_dag()
        dp = "data_processing"
        # ingest_load depends on ingest_raw
        assert (dp, "ingest_raw") in dag[(dp, "ingest_load")].depends_on
        # qc_filter depends on ingest_load
        assert (dp, "ingest_load") in dag[(dp, "qc_filter")].depends_on
        # metadata_attach depends on qc_filter
        assert (dp, "qc_filter") in dag[(dp, "metadata_attach")].depends_on
        # preprocess depends on metadata_attach
        assert (dp, "metadata_attach") in dag[(dp, "preprocess")].depends_on
        # demographics depends on preprocess
        assert (dp, "preprocess") in dag[(dp, "demographics")].depends_on
        # scoring depends on preprocess
        assert (dp, "preprocess") in dag[(dp, "scoring")].depends_on
        # celltype_manual depends on scoring
        assert (dp, "scoring") in dag[(dp, "celltype_manual")].depends_on
        # biology depends on scoring
        assert (dp, "scoring") in dag[(dp, "biology")].depends_on
        # meta_analysis depends on biology
        assert (dp, "biology") in dag[(dp, "meta_analysis")].depends_on


# ---------------------------------------------------------------------------
# Tuple-based get_available_next
# ---------------------------------------------------------------------------


class TestTupleGetAvailableNext:
    def test_empty_completed_returns_root_tuples(self):
        from sc_tools.pipeline import get_available_next

        available = get_available_next([])
        assert ("data_processing", "ingest_raw") in available
        assert ("data_processing", "qc_filter") not in available

    def test_after_ingest_raw_tuple(self):
        from sc_tools.pipeline import get_available_next

        available = get_available_next([("data_processing", "ingest_raw")])
        assert ("data_processing", "ingest_load") in available

    def test_branching_at_preprocess_tuples(self):
        from sc_tools.pipeline import get_available_next

        dp = "data_processing"
        completed = [
            (dp, "ingest_raw"),
            (dp, "ingest_load"),
            (dp, "qc_filter"),
            (dp, "metadata_attach"),
            (dp, "preprocess"),
        ]
        available = get_available_next(completed)
        assert (dp, "demographics") in available
        assert (dp, "scoring") in available

    def test_iterative_phase_re_available_tuples(self):
        from sc_tools.pipeline import get_available_next

        dp = "data_processing"
        completed = [
            (dp, "ingest_raw"),
            (dp, "ingest_load"),
            (dp, "qc_filter"),
            (dp, "metadata_attach"),
            (dp, "preprocess"),
            (dp, "scoring"),
            (dp, "celltype_manual"),
        ]
        available = get_available_next(completed)
        # celltype_manual is iterative -> still available
        assert (dp, "celltype_manual") in available

    def test_all_phases_complete_tuples(self):
        from sc_tools.pipeline import get_available_next, get_dag

        dp = "data_processing"
        all_phases = [
            (dp, "ingest_raw"),
            (dp, "ingest_load"),
            (dp, "qc_filter"),
            (dp, "metadata_attach"),
            (dp, "preprocess"),
            (dp, "demographics"),
            (dp, "scoring"),
            (dp, "celltype_manual"),
            (dp, "biology"),
            (dp, "meta_analysis"),
        ]
        available = get_available_next(all_phases)
        dag = get_dag()
        for p in available:
            assert dag[p].iterative, f"Non-iterative phase {p!r} appeared after full completion"


# ---------------------------------------------------------------------------
# flat_slug_to_tuple / tuple_to_display helpers
# ---------------------------------------------------------------------------


class TestHelpers:
    def test_flat_slug_to_tuple_standard(self):
        from sc_tools.pipeline import flat_slug_to_tuple

        assert flat_slug_to_tuple("qc_filter") == ("data_processing", "qc_filter")
        assert flat_slug_to_tuple("ingest_raw") == ("data_processing", "ingest_raw")
        assert flat_slug_to_tuple("meta_analysis") == ("data_processing", "meta_analysis")

    def test_flat_slug_to_tuple_unknown_raises(self):
        from sc_tools.pipeline import flat_slug_to_tuple

        with pytest.raises(KeyError, match="not registered"):
            flat_slug_to_tuple("nonexistent_slug")

    def test_tuple_to_display(self):
        from sc_tools.pipeline import tuple_to_display

        assert tuple_to_display(("data_processing", "qc_filter")) == "data_processing/qc_filter"
        assert tuple_to_display(("discovery", "clustering_v1")) == "discovery/clustering_v1"

    def test_flat_slug_to_tuple_all_standard(self):
        """All 10 standard slugs should map correctly."""
        from sc_tools.pipeline import STANDARD_PHASES, flat_slug_to_tuple

        for slug in STANDARD_PHASES:
            result = flat_slug_to_tuple(slug)
            assert result == ("data_processing", slug)


# ---------------------------------------------------------------------------
# Backward compat: flat slugs still work with get_available_next
# ---------------------------------------------------------------------------


class TestBackwardCompatFlatSlugs:
    def test_get_available_next_accepts_flat_slugs(self):
        """Old code passing flat slugs should still work."""
        from sc_tools.pipeline import get_available_next

        available = get_available_next(["ingest_raw"])
        # Should return tuples even when given flat slugs
        assert ("data_processing", "ingest_load") in available

    def test_get_available_next_mixed_input(self):
        """Accept both flat slugs and tuples in the same list."""
        from sc_tools.pipeline import get_available_next

        completed = [
            "ingest_raw",
            ("data_processing", "ingest_load"),
        ]
        available = get_available_next(completed)
        assert ("data_processing", "qc_filter") in available

    def test_get_phase_accepts_flat_slug(self):
        """get_phase should still accept flat slugs for backward compat."""
        from sc_tools.pipeline import get_phase

        spec = get_phase("qc_filter")
        assert spec.label == "QC Filtering + Concatenation"

    def test_get_phase_accepts_tuple(self):
        """get_phase should also accept tuples."""
        from sc_tools.pipeline import get_phase

        spec = get_phase(("data_processing", "qc_filter"))
        assert spec.label == "QC Filtering + Concatenation"

    def test_get_phase_checkpoint_accepts_flat_slug(self):
        """get_phase_checkpoint should still accept flat slugs."""
        from sc_tools.pipeline import get_phase_checkpoint

        cp = get_phase_checkpoint("qc_filter")
        assert cp == "results/adata.filtered.h5ad"

    def test_get_phase_checkpoint_accepts_tuple(self):
        """get_phase_checkpoint should also accept tuples."""
        from sc_tools.pipeline import get_phase_checkpoint

        cp = get_phase_checkpoint(("data_processing", "qc_filter"))
        assert cp == "results/adata.filtered.h5ad"
