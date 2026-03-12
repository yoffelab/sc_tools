"""Unit tests for sc_tools.registry."""

from __future__ import annotations

import json
import warnings

import pytest

# Skip all tests if SQLAlchemy is not installed
sqlalchemy = pytest.importorskip("sqlalchemy", reason="sqlalchemy not installed")

# Suppress register_dataset deprecation warnings in existing tests
pytestmark = pytest.mark.filterwarnings("ignore::DeprecationWarning")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def reg(tmp_path):
    """Return an in-memory SQLite Registry for testing."""
    from sc_tools.registry import Registry

    db_url = f"sqlite:///{tmp_path / 'test_registry.db'}"
    return Registry(db_url=db_url)


# ---------------------------------------------------------------------------
# Projects
# ---------------------------------------------------------------------------


class TestProjects:
    def test_add_project(self, reg):
        pid = reg.add_project("proj_a", platform="visium", data_type="visium")
        assert isinstance(pid, int)
        assert pid > 0

    def test_add_project_idempotent(self, reg):
        pid1 = reg.add_project("proj_a", platform="visium")
        pid2 = reg.add_project("proj_a", platform="visium")
        assert pid1 == pid2

    def test_get_project(self, reg):
        reg.add_project("proj_b", platform="imc")
        proj = reg.get_project("proj_b")
        assert proj is not None
        assert proj["name"] == "proj_b"
        assert proj["platform"] == "imc"
        assert proj["status"] == "active"

    def test_get_project_missing(self, reg):
        assert reg.get_project("nonexistent") is None

    def test_list_projects(self, reg):
        reg.add_project("p1", platform="visium")
        reg.add_project("p2", platform="xenium")
        projects = reg.list_projects()
        names = [p["name"] for p in projects]
        assert "p1" in names
        assert "p2" in names

    def test_mark_phase_complete(self, reg):
        reg.add_project("proj_c", platform="visium")
        reg.mark_phase_complete("proj_c", "p1")
        reg.mark_phase_complete("proj_c", "p2")
        proj = reg.get_project("proj_c")
        phases = json.loads(proj["phases_complete"])
        assert "p1" in phases
        assert "p2" in phases

    def test_mark_phase_complete_idempotent(self, reg):
        reg.add_project("proj_d", platform="visium")
        reg.mark_phase_complete("proj_d", "p1")
        reg.mark_phase_complete("proj_d", "p1")  # duplicate
        proj = reg.get_project("proj_d")
        phases = json.loads(proj["phases_complete"])
        assert phases.count("p1") == 1

    def test_mark_phase_complete_missing_project(self, reg):
        with pytest.raises(ValueError, match="not found"):
            reg.mark_phase_complete("nonexistent", "p1")


# ---------------------------------------------------------------------------
# Datasets
# ---------------------------------------------------------------------------


class TestDatasets:
    def test_register_dataset(self, reg):
        reg.add_project("proj", platform="visium")
        ds_id = reg.register_dataset("proj", phase="qc_filter", uri="/data/adata.raw.h5ad")
        assert isinstance(ds_id, int)

    def test_register_dataset_missing_project(self, reg):
        with pytest.raises(ValueError, match="not found"):
            reg.register_dataset("nonexistent", phase="p1", uri="/x.h5ad")

    def test_get_dataset_uri(self, reg):
        reg.add_project("proj", platform="visium")
        reg.register_dataset("proj", phase="preprocess", uri="/results/adata.normalized.h5ad")
        uri = reg.get_dataset_uri("proj", phase="preprocess")
        assert uri == "/results/adata.normalized.h5ad"

    def test_get_dataset_uri_with_sample(self, reg):
        reg.add_project("proj", platform="imc")
        reg.register_dataset("proj", phase="p0", uri="/data/s1/adata.p0.h5ad", sample_id="s1")
        reg.register_dataset("proj", phase="p0", uri="/data/s2/adata.p0.h5ad", sample_id="s2")
        assert reg.get_dataset_uri("proj", "p0", sample_id="s1") == "/data/s1/adata.p0.h5ad"
        assert reg.get_dataset_uri("proj", "p0", sample_id="s2") == "/data/s2/adata.p0.h5ad"

    def test_get_dataset_uri_missing(self, reg):
        reg.add_project("proj", platform="visium")
        assert reg.get_dataset_uri("proj", "p99") is None

    def test_get_dataset_uri_missing_project(self, reg):
        assert reg.get_dataset_uri("nonexistent", "p1") is None

    def test_list_datasets(self, reg):
        reg.add_project("proj", platform="visium")
        reg.register_dataset("proj", phase="p1", uri="/a.h5ad")
        reg.register_dataset("proj", phase="p3", uri="/b.h5ad")
        datasets = reg.list_datasets(project_name="proj")
        assert len(datasets) == 2
        phases = {d["phase"] for d in datasets}
        assert phases == {"p1", "p3"}

    def test_list_datasets_phase_filter(self, reg):
        reg.add_project("proj", platform="visium")
        reg.register_dataset("proj", phase="p1", uri="/a.h5ad")
        reg.register_dataset("proj", phase="p3", uri="/b.h5ad")
        datasets = reg.list_datasets(project_name="proj", phase="p1")
        assert len(datasets) == 1
        assert datasets[0]["phase"] == "p1"

    def test_list_datasets_all(self, reg):
        reg.add_project("p1", platform="visium")
        reg.add_project("p2", platform="imc")
        reg.register_dataset("p1", phase="p1", uri="/a.h5ad")
        reg.register_dataset("p2", phase="p1", uri="/b.h5ad")
        all_ds = reg.list_datasets()
        assert len(all_ds) == 2

    def test_update_dataset_status(self, reg):
        reg.add_project("proj", platform="visium")
        ds_id = reg.register_dataset("proj", phase="p1", uri="/a.h5ad", status="pending")
        reg.update_dataset_status(ds_id, "ready")
        datasets = reg.list_datasets(project_name="proj")
        assert datasets[0]["status"] == "ready"

    def test_update_dataset_status_missing(self, reg):
        with pytest.raises(ValueError, match="not found"):
            reg.update_dataset_status(9999, "ready")


# ---------------------------------------------------------------------------
# SLURM jobs
# ---------------------------------------------------------------------------


class TestSlurmJobs:
    def test_register_job(self, reg):
        reg.add_project("proj", platform="visium")
        job_id = reg.register_job(
            "proj", sample_id="s1", phase="p0", cluster="brb", slurm_job_id="12345"
        )
        assert isinstance(job_id, int)

    def test_register_job_missing_project(self, reg):
        with pytest.raises(ValueError, match="not found"):
            reg.register_job("nonexistent", None, "p0", "brb", "123")

    def test_update_job_status_completed(self, reg):
        reg.add_project("proj", platform="visium")
        reg.register_job("proj", None, "p0", "brb", "999")
        reg.update_job_status("999", "completed")
        jobs = reg.list_active_jobs()
        assert not any(j["slurm_job_id"] == "999" for j in jobs)

    def test_update_job_status_with_error(self, reg):
        reg.add_project("proj", platform="visium")
        reg.register_job("proj", None, "p0", "cayuga", "888")
        reg.update_job_status("888", "failed", error="OOM")
        # Still findable — it's not active
        jobs = reg.list_active_jobs()
        assert not any(j["slurm_job_id"] == "888" for j in jobs)

    def test_update_job_missing(self, reg):
        with pytest.raises(ValueError, match="not found"):
            reg.update_job_status("nonexistent_id", "completed")

    def test_list_active_jobs(self, reg):
        reg.add_project("proj", platform="visium")
        reg.register_job("proj", None, "p0", "brb", "111")
        reg.register_job("proj", None, "p0", "brb", "222")
        active = reg.list_active_jobs()
        slurm_ids = {j["slurm_job_id"] for j in active}
        assert "111" in slurm_ids
        assert "222" in slurm_ids

        reg.update_job_status("111", "completed")
        active = reg.list_active_jobs()
        assert not any(j["slurm_job_id"] == "111" for j in active)


# ---------------------------------------------------------------------------
# Agent tasks
# ---------------------------------------------------------------------------


class TestAgentTasks:
    def test_start_and_finish_task(self, reg):
        reg.add_project("proj", platform="visium")
        task_id = reg.start_task("qc_report", "proj", inputs={"adata": "/a.h5ad"})
        assert isinstance(task_id, int)

        running = reg.list_running_tasks()
        assert any(t["id"] == task_id for t in running)

        reg.finish_task(task_id, outputs={"report": "/qc.html"})
        running = reg.list_running_tasks()
        assert not any(t["id"] == task_id for t in running)

    def test_finish_task_with_error(self, reg):
        reg.add_project("proj", platform="visium")
        task_id = reg.start_task("scoring", "proj")
        reg.finish_task(task_id, error="segfault")
        running = reg.list_running_tasks()
        assert not any(t["id"] == task_id for t in running)

    def test_start_task_missing_project(self, reg):
        with pytest.raises(ValueError, match="not found"):
            reg.start_task("qc", "nonexistent")

    def test_finish_task_missing(self, reg):
        with pytest.raises(ValueError, match="not found"):
            reg.finish_task(9999)


# ---------------------------------------------------------------------------
# Status summary
# ---------------------------------------------------------------------------


class TestStatus:
    def test_status_empty(self, reg):
        s = reg.status()
        assert s["n_projects"] == 0
        assert s["n_datasets"] == 0
        assert s["active_slurm_jobs"] == 0
        assert s["running_agent_tasks"] == 0

    def test_status_with_data(self, reg):
        reg.add_project("proj", platform="visium")
        reg.register_dataset("proj", phase="p1", uri="/a.h5ad")
        reg.register_job("proj", None, "p0", "brb", "42")
        reg.start_task("qc", "proj")
        s = reg.status()
        assert s["n_projects"] == 1
        assert s["n_datasets"] == 1
        assert s["active_slurm_jobs"] == 1
        assert s["running_agent_tasks"] == 1
        assert "proj" in s["active_projects"]


# ---------------------------------------------------------------------------
# Technology taxonomy
# ---------------------------------------------------------------------------


class TestProjectTaxonomy:
    def test_add_project_with_domain(self, reg):
        pid = reg.add_project(
            "ggo_visium",
            platform="visium",
            data_type="visium",
            domain="spatial_transcriptomics",
            imaging_modality="sequencing_based",
        )
        proj = reg.get_project("ggo_visium")
        assert proj is not None
        assert proj["domain"] == "spatial_transcriptomics"
        assert proj["imaging_modality"] == "sequencing_based"
        assert pid > 0

    def test_add_project_imc_taxonomy(self, reg):
        reg.add_project(
            "lymph_dlbcl",
            platform="imc",
            domain="spatial_proteomics",
            imaging_modality="mass_spec_imaging",
        )
        proj = reg.get_project("lymph_dlbcl")
        assert proj["domain"] == "spatial_proteomics"
        assert proj["imaging_modality"] == "mass_spec_imaging"

    def test_add_project_null_domain_ok(self, reg):
        reg.add_project("no_domain_proj", platform="visium")
        proj = reg.get_project("no_domain_proj")
        assert proj is not None
        assert proj["domain"] is None
        assert proj["imaging_modality"] is None

    def test_add_project_idempotent_with_taxonomy(self, reg):
        pid1 = reg.add_project("proj", platform="imc", domain="spatial_proteomics")
        pid2 = reg.add_project("proj", platform="imc", domain="spatial_proteomics")
        assert pid1 == pid2

    def test_list_projects_includes_taxonomy(self, reg):
        reg.add_project("p1", platform="visium", domain="spatial_transcriptomics")
        reg.add_project("p2", platform="imc", domain="spatial_proteomics")
        projects = reg.list_projects()
        domains = {p["name"]: p["domain"] for p in projects}
        assert domains["p1"] == "spatial_transcriptomics"
        assert domains["p2"] == "spatial_proteomics"


# ---------------------------------------------------------------------------
# Dataset file role, validated, n_obs/n_vars
# ---------------------------------------------------------------------------


class TestDatasetFileRole:
    def test_register_primary_dataset(self, reg):
        reg.add_project("proj", platform="visium")
        reg.register_dataset(
            "proj",
            phase="qc_filter",
            uri="/results/adata.raw.h5ad",
            file_role="primary",
            n_obs=50000,
            n_vars=3000,
        )
        datasets = reg.list_datasets("proj")
        assert len(datasets) == 1
        assert datasets[0]["file_role"] == "primary"
        assert datasets[0]["n_obs"] == 50000
        assert datasets[0]["n_vars"] == 3000

    def test_register_supplementary_spatialdata(self, reg):
        reg.add_project("proj", platform="visium_hd")
        reg.register_dataset(
            "proj",
            phase="ingest_load",
            uri="/data/s1/adata.h5ad",
            file_role="primary",
            sample_id="s1",
        )
        reg.register_dataset(
            "proj",
            phase="ingest_load",
            uri="/data/s1/spatialdata.zarr",
            file_role="spatialdata",
            sample_id="s1",
            fmt="zarr",
        )
        datasets = reg.list_datasets("proj")
        assert len(datasets) == 2
        roles = {d["file_role"] for d in datasets}
        assert roles == {"primary", "spatialdata"}

    def test_entry_point_role(self, reg):
        reg.add_project("proj", platform="imc")
        reg.register_dataset(
            "proj",
            phase="celltype_manual",
            uri="/results/adata.celltyped.h5ad",
            file_role="entry_point",
        )
        ds = reg.list_datasets("proj")[0]
        assert ds["file_role"] == "entry_point"

    def test_validated_flag_default_false(self, reg):
        reg.add_project("proj", platform="visium")
        reg.register_dataset("proj", phase="qc_filter", uri="/a.h5ad")
        ds = reg.list_datasets("proj")[0]
        assert ds["validated"] is False

    def test_mark_dataset_validated(self, reg):
        reg.add_project("proj", platform="visium")
        ds_id = reg.register_dataset("proj", phase="qc_filter", uri="/a.h5ad")
        reg.mark_dataset_validated(ds_id)
        ds = reg.list_datasets("proj")[0]
        assert ds["validated"] is True

    def test_mark_dataset_validated_missing(self, reg):
        with pytest.raises(ValueError, match="not found"):
            reg.mark_dataset_validated(9999)

    def test_register_dataset_with_all_new_fields(self, reg):
        reg.add_project("proj", platform="imc")
        reg.register_dataset(
            "proj",
            phase="scoring",
            uri="/results/adata.scored.h5ad",
            file_role="primary",
            validated=True,
            n_obs=248000,
            n_vars=38,
        )
        ds = reg.list_datasets("proj")[0]
        assert ds["file_role"] == "primary"
        assert ds["validated"] is True
        assert ds["n_obs"] == 248000
        assert ds["n_vars"] == 38


# ---------------------------------------------------------------------------
# Project phases (upsert_phase, get_phase, list_phases)
# ---------------------------------------------------------------------------


class TestProjectPhases:
    def test_upsert_phase_creates_new(self, reg):
        reg.add_project("proj", platform="visium")
        reg.upsert_phase("proj", "qc_filter", status="in_progress", n_obs=50000, n_samples=8)
        row = reg.get_phase("proj", "qc_filter")
        assert row is not None
        assert row["status"] == "in_progress"
        assert row["n_obs"] == 50000
        assert row["n_samples"] == 8

    def test_upsert_phase_updates_existing(self, reg):
        reg.add_project("proj", platform="visium")
        reg.upsert_phase("proj", "qc_filter", status="in_progress")
        reg.upsert_phase("proj", "qc_filter", status="complete", n_obs=50000)
        row = reg.get_phase("proj", "qc_filter")
        assert row["status"] == "complete"
        assert row["n_obs"] == 50000

    def test_upsert_phase_preserves_existing_n_obs_if_not_provided(self, reg):
        reg.add_project("proj", platform="visium")
        reg.upsert_phase("proj", "qc_filter", status="in_progress", n_obs=50000)
        reg.upsert_phase("proj", "qc_filter", status="complete")
        row = reg.get_phase("proj", "qc_filter")
        assert row["n_obs"] == 50000  # not overwritten by None

    def test_get_phase_missing_returns_none(self, reg):
        reg.add_project("proj", platform="visium")
        assert reg.get_phase("proj", "nonexistent") is None

    def test_get_phase_missing_project_returns_none(self, reg):
        assert reg.get_phase("nonexistent", "qc_filter") is None

    def test_list_phases_empty(self, reg):
        reg.add_project("proj", platform="visium")
        assert reg.list_phases("proj") == []

    def test_list_phases_multiple(self, reg):
        reg.add_project("proj", platform="imc")
        reg.upsert_phase("proj", "qc_filter", status="complete")
        reg.upsert_phase("proj", "preprocess", status="complete")
        reg.upsert_phase("proj", "scoring", status="in_progress")
        phases = reg.list_phases("proj")
        assert len(phases) == 3
        slugs = [p["phase"] for p in phases]
        # ordered by phase slug (alphabetical)
        assert "preprocess" in slugs
        assert "qc_filter" in slugs
        assert "scoring" in slugs

    def test_list_phases_missing_project(self, reg):
        assert reg.list_phases("nonexistent") == []

    def test_mark_phase_complete_syncs_project_phases(self, reg):
        reg.add_project("proj", platform="visium")
        reg.mark_phase_complete("proj", "qc_filter")
        row = reg.get_phase("proj", "qc_filter")
        assert row is not None
        assert row["status"] == "complete"
        # Legacy phases_complete should also be updated
        proj = reg.get_project("proj")
        phases = json.loads(proj["phases_complete"])
        assert "qc_filter" in phases

    def test_entry_phase_flag(self, reg):
        reg.add_project("proj", platform="imc")
        reg.upsert_phase("proj", "celltype_manual", status="in_progress", entry_phase=True)
        row = reg.get_phase("proj", "celltype_manual")
        assert row["entry_phase"] is True

    def test_upsert_phase_with_primary_dataset_id(self, reg):
        reg.add_project("proj", platform="visium")
        ds_id = reg.register_dataset("proj", phase="qc_filter", uri="/a.h5ad")
        reg.upsert_phase("proj", "qc_filter", status="complete", primary_dataset_id=ds_id)
        row = reg.get_phase("proj", "qc_filter")
        assert row["primary_dataset_id"] == ds_id

    def test_upsert_phase_notes(self, reg):
        reg.add_project("proj", platform="imc")
        reg.upsert_phase(
            "proj",
            "preprocess",
            status="complete",
            notes="9-method benchmark; Z-score+Harmony selected",
        )
        row = reg.get_phase("proj", "preprocess")
        assert "Z-score+Harmony" in row["notes"]

    def test_upsert_phase_missing_project_raises(self, reg):
        with pytest.raises(ValueError, match="not found"):
            reg.upsert_phase("nonexistent", "qc_filter")

    def test_status_includes_phase_summary(self, reg):
        reg.add_project("proj", platform="visium")
        reg.upsert_phase("proj", "qc_filter", status="complete")
        reg.upsert_phase("proj", "preprocess", status="in_progress")
        s = reg.status()
        assert "proj" in s["active_projects"]
        summary = s.get("phase_summary", {}).get("proj", {})
        assert summary.get("complete", 0) == 1
        assert summary.get("in_progress", 0) == 1


# ---------------------------------------------------------------------------
# Subjects
# ---------------------------------------------------------------------------


class TestSubjects:
    def test_add_subject(self, reg):
        sid = reg.add_subject("PT001", organism="human", sex="M", diagnosis="DLBCL")
        assert isinstance(sid, int)
        assert sid > 0

    def test_add_subject_idempotent(self, reg):
        sid1 = reg.add_subject("PT001")
        sid2 = reg.add_subject("PT001")
        assert sid1 == sid2

    def test_get_subject(self, reg):
        reg.add_subject("PT002", organism="mouse", sex="F", diagnosis="GGO")
        subj = reg.get_subject("PT002")
        assert subj is not None
        assert subj["subject_id"] == "PT002"
        assert subj["organism"] == "mouse"
        assert subj["sex"] == "F"
        assert subj["diagnosis"] == "GGO"

    def test_get_subject_missing(self, reg):
        assert reg.get_subject("NONEXISTENT") is None

    def test_list_subjects_all(self, reg):
        reg.add_subject("PT001")
        reg.add_subject("PT002")
        subjects = reg.list_subjects()
        ids = [s["subject_id"] for s in subjects]
        assert "PT001" in ids
        assert "PT002" in ids

    def test_list_subjects_filter_diagnosis(self, reg):
        reg.add_subject("PT001", diagnosis="DLBCL")
        reg.add_subject("PT002", diagnosis="UC")
        result = reg.list_subjects(diagnosis="DLBCL")
        assert len(result) == 1
        assert result[0]["subject_id"] == "PT001"

    def test_list_subjects_filter_tissue(self, reg):
        reg.add_subject("PT001", tissue_of_origin="lymph_node")
        reg.add_subject("PT002", tissue_of_origin="colon")
        result = reg.list_subjects(tissue="lymph")
        assert len(result) == 1
        assert result[0]["subject_id"] == "PT001"

    def test_list_subjects_by_project(self, reg):
        reg.add_project("proj_a", platform="visium")
        reg.add_subject("PT001")
        reg.add_subject("PT002")
        reg.link_subject_to_project("PT001", "proj_a")
        result = reg.list_subjects(project_name="proj_a")
        assert len(result) == 1
        assert result[0]["subject_id"] == "PT001"

    def test_add_subject_with_clinical_overflow(self, reg):
        reg.add_subject(
            "PT003",
            diagnosis="DLBCL",
            clinical_metadata_json='{"response": "CR", "cycles": 6}',
        )
        subj = reg.get_subject("PT003")
        assert subj["clinical_metadata_json"] is not None
        import json

        data = json.loads(subj["clinical_metadata_json"])
        assert data["response"] == "CR"

    def test_link_subject_to_project(self, reg):
        reg.add_project("proj_a", platform="visium")
        reg.add_subject("PT001")
        link_id = reg.link_subject_to_project("PT001", "proj_a")
        assert isinstance(link_id, int)

    def test_link_subject_idempotent(self, reg):
        reg.add_project("proj_a", platform="visium")
        reg.add_subject("PT001")
        id1 = reg.link_subject_to_project("PT001", "proj_a")
        id2 = reg.link_subject_to_project("PT001", "proj_a")
        assert id1 == id2

    def test_link_subject_missing_subject_raises(self, reg):
        reg.add_project("proj_a", platform="visium")
        with pytest.raises(ValueError, match="not found"):
            reg.link_subject_to_project("NONEXISTENT", "proj_a")

    def test_link_subject_missing_project_raises(self, reg):
        reg.add_subject("PT001")
        with pytest.raises(ValueError, match="not found"):
            reg.link_subject_to_project("PT001", "nonexistent_proj")


# ---------------------------------------------------------------------------
# Samples
# ---------------------------------------------------------------------------


class TestSamples:
    def test_add_sample(self, reg):
        reg.add_project("proj", platform="visium")
        reg.add_subject("PT001")
        sid = reg.add_sample("S001", subject_id="PT001", project_name="proj", tissue="colon")
        assert isinstance(sid, int)
        assert sid > 0

    def test_add_sample_without_subject(self, reg):
        reg.add_project("proj", platform="visium")
        sid = reg.add_sample("S002", project_name="proj")
        assert isinstance(sid, int)

    def test_add_sample_missing_subject_raises(self, reg):
        reg.add_project("proj", platform="visium")
        with pytest.raises(ValueError, match="not found"):
            reg.add_sample("S001", subject_id="NONEXISTENT", project_name="proj")

    def test_add_sample_missing_project_raises(self, reg):
        with pytest.raises(ValueError, match="not found"):
            reg.add_sample("S001", project_name="nonexistent")

    def test_get_sample(self, reg):
        reg.add_project("proj", platform="imc")
        reg.add_subject("PT001")
        reg.add_sample(
            "S001",
            subject_id="PT001",
            project_name="proj",
            tissue="lymph_node",
            fixation_method="FFPE",
        )
        sample = reg.get_sample("S001")
        assert sample is not None
        assert sample["tissue"] == "lymph_node"
        assert sample["fixation_method"] == "FFPE"

    def test_get_sample_scoped_to_project(self, reg):
        reg.add_project("proj_a", platform="visium")
        reg.add_project("proj_b", platform="imc")
        reg.add_sample("S001", project_name="proj_a")
        reg.add_sample("S001", project_name="proj_b")
        sample_a = reg.get_sample("S001", project_name="proj_a")
        sample_b = reg.get_sample("S001", project_name="proj_b")
        assert sample_a["project_id"] != sample_b["project_id"]

    def test_get_sample_missing(self, reg):
        assert reg.get_sample("NONEXISTENT") is None

    def test_list_samples_by_project(self, reg):
        reg.add_project("proj", platform="visium")
        reg.add_sample("S001", project_name="proj")
        reg.add_sample("S002", project_name="proj")
        samples = reg.list_samples(project_name="proj")
        assert len(samples) == 2

    def test_list_samples_by_subject(self, reg):
        reg.add_project("proj", platform="visium")
        reg.add_subject("PT001")
        reg.add_subject("PT002")
        reg.add_sample("S001", subject_id="PT001", project_name="proj")
        reg.add_sample("S002", subject_id="PT002", project_name="proj")
        samples = reg.list_samples(subject_id="PT001")
        assert len(samples) == 1
        assert samples[0]["sample_id"] == "S001"

    def test_list_samples_by_batch(self, reg):
        reg.add_project("proj", platform="visium")
        reg.add_sample("S001", project_name="proj", batch="batch1")
        reg.add_sample("S002", project_name="proj", batch="batch2")
        result = reg.list_samples(batch="batch1")
        assert len(result) == 1


# ---------------------------------------------------------------------------
# BioData CRUD and polymorphism
# ---------------------------------------------------------------------------


class TestBioData:
    def test_register_spatial_seq(self, reg):
        reg.add_project("proj", platform="visium")
        bd_id = reg.register_biodata(
            "proj",
            "spatial_seq",
            "visium",
            "/results/adata.filtered.h5ad",
            fmt="h5ad",
            phase="qc_filter",
            n_obs=50000,
        )
        assert isinstance(bd_id, int)

    def test_register_image(self, reg):
        reg.add_project("proj", platform="imc")
        bd_id = reg.register_biodata(
            "proj",
            "image",
            "imc",
            "/data/sample1.tiff",
            fmt="tiff",
            image_type="multiplexed",
            n_channels=40,
        )
        assert isinstance(bd_id, int)

    def test_register_rnaseq(self, reg):
        reg.add_project("proj", platform="chromium_3p")
        bd_id = reg.register_biodata(
            "proj",
            "rnaseq",
            "chromium_3p",
            "/data/counts.h5ad",
            fmt="h5ad",
            chemistry="chromium_v3",
        )
        assert isinstance(bd_id, int)

    def test_register_epigenomics(self, reg):
        reg.add_project("proj", platform="atac_seq")
        bd_id = reg.register_biodata(
            "proj",
            "epigenomics",
            "atac_seq",
            "/data/peaks.bed",
            fmt="bed",
            assay_type="atac_seq",
            n_peaks=100000,
        )
        assert isinstance(bd_id, int)

    def test_register_genome_seq(self, reg):
        reg.add_project("proj", platform="illumina_wgs")
        bd_id = reg.register_biodata(
            "proj",
            "genome_seq",
            "illumina_wgs",
            "/data/aligned.bam",
            fmt="bam",
            sequencing_type="wgs",
            coverage_mean=30.0,
        )
        assert isinstance(bd_id, int)

    def test_register_missing_project_raises(self, reg):
        with pytest.raises(ValueError, match="not found"):
            reg.register_biodata("nonexistent", "spatial_seq", "visium", "/a.h5ad")

    def test_register_invalid_category_raises(self, reg):
        reg.add_project("proj", platform="visium")
        with pytest.raises(ValueError, match="Unknown BioData category"):
            reg.register_biodata("proj", "invalid_cat", "visium", "/a.h5ad")

    def test_get_biodata(self, reg):
        reg.add_project("proj", platform="visium")
        bd_id = reg.register_biodata(
            "proj",
            "spatial_seq",
            "visium",
            "/results/adata.h5ad",
            fmt="h5ad",
            spatial_resolution="spot",
            panel_size=None,
        )
        bd = reg.get_biodata(bd_id)
        assert bd is not None
        assert bd["type"] == "spatial_seq"
        assert bd["platform"] == "visium"
        assert bd["uri"] == "/results/adata.h5ad"
        # Check child columns are present
        assert "spatial_resolution" in bd

    def test_get_biodata_missing(self, reg):
        assert reg.get_biodata(9999) is None

    def test_list_biodata(self, reg):
        reg.add_project("proj", platform="visium")
        reg.register_biodata("proj", "spatial_seq", "visium", "/a.h5ad")
        reg.register_biodata("proj", "spatial_seq", "xenium", "/b.h5ad")
        result = reg.list_biodata(project_name="proj")
        assert len(result) == 2

    def test_list_biodata_filter_category(self, reg):
        reg.add_project("proj", platform="visium")
        reg.register_biodata("proj", "spatial_seq", "visium", "/a.h5ad")
        reg.register_biodata("proj", "image", "imc", "/b.tiff", fmt="tiff")
        result = reg.list_biodata(project_name="proj", category="image")
        assert len(result) == 1
        assert result[0]["platform"] == "imc"

    def test_list_biodata_filter_platform(self, reg):
        reg.add_project("proj", platform="visium")
        reg.register_biodata("proj", "spatial_seq", "visium", "/a.h5ad")
        reg.register_biodata("proj", "spatial_seq", "xenium", "/b.h5ad")
        result = reg.list_biodata(project_name="proj", platform="xenium")
        assert len(result) == 1

    def test_list_biodata_filter_phase(self, reg):
        reg.add_project("proj", platform="visium")
        reg.register_biodata("proj", "spatial_seq", "visium", "/a.h5ad", phase="qc_filter")
        reg.register_biodata("proj", "spatial_seq", "visium", "/b.h5ad", phase="preprocess")
        result = reg.list_biodata(project_name="proj", phase="qc_filter")
        assert len(result) == 1


class TestBioDataPolymorphism:
    def test_image_type_specific_fields(self, reg):
        reg.add_project("proj", platform="imc")
        bd_id = reg.register_biodata(
            "proj",
            "image",
            "imc",
            "/data/tiff.tiff",
            fmt="tiff",
            image_type="multiplexed",
            n_channels=40,
            pixel_size_um=1.0,
            staining_protocol="IMC",
        )
        bd = reg.get_biodata(bd_id)
        assert bd["image_type"] == "multiplexed"
        assert bd["n_channels"] == 40
        assert bd["pixel_size_um"] == 1.0
        assert bd["staining_protocol"] == "IMC"

    def test_spatial_seq_type_specific_fields(self, reg):
        reg.add_project("proj", platform="visium_hd")
        bd_id = reg.register_biodata(
            "proj",
            "spatial_seq",
            "visium_hd",
            "/data/bins.h5ad",
            fmt="h5ad",
            spatial_resolution="spot",
            bin_size_um=8.0,
        )
        bd = reg.get_biodata(bd_id)
        assert bd["spatial_resolution"] == "spot"
        assert bd["bin_size_um"] == 8.0

    def test_rnaseq_type_specific_fields(self, reg):
        reg.add_project("proj", platform="chromium_3p")
        bd_id = reg.register_biodata(
            "proj",
            "rnaseq",
            "chromium_3p",
            "/data/counts.h5ad",
            fmt="h5ad",
            chemistry="chromium_v3",
            reference_genome="GRCh38",
        )
        bd = reg.get_biodata(bd_id)
        assert bd["chemistry"] == "chromium_v3"
        assert bd["reference_genome"] == "GRCh38"

    def test_epigenomics_type_specific_fields(self, reg):
        reg.add_project("proj", platform="chip_seq")
        bd_id = reg.register_biodata(
            "proj",
            "epigenomics",
            "chip_seq",
            "/data/peaks.bed",
            fmt="bed",
            assay_type="chip_seq",
            target_protein="H3K27ac",
            n_peaks=50000,
        )
        bd = reg.get_biodata(bd_id)
        assert bd["assay_type"] == "chip_seq"
        assert bd["target_protein"] == "H3K27ac"
        assert bd["n_peaks"] == 50000

    def test_genome_seq_type_specific_fields(self, reg):
        reg.add_project("proj", platform="pacbio_hifi")
        bd_id = reg.register_biodata(
            "proj",
            "genome_seq",
            "pacbio_hifi",
            "/data/aligned.bam",
            fmt="bam",
            sequencing_type="wgs",
            coverage_mean=30.0,
            reference_genome="GRCh38",
        )
        bd = reg.get_biodata(bd_id)
        assert bd["sequencing_type"] == "wgs"
        assert bd["coverage_mean"] == 30.0

    def test_auto_fill_from_platform_registry(self, reg):
        reg.add_project("proj", platform="xenium")
        bd_id = reg.register_biodata(
            "proj",
            "spatial_seq",
            "xenium",
            "/data/xenium.h5ad",
            fmt="h5ad",
        )
        bd = reg.get_biodata(bd_id)
        # Should have been auto-filled from KNOWN_PLATFORMS
        assert bd["measurement"] == "rna"
        assert bd["spatial"] is True
        assert bd["resolution"] == "subcellular"


class TestProjectDataSummary:
    def test_empty_project(self, reg):
        reg.add_project("proj", platform="visium")
        summary = reg.project_data_summary("proj")
        assert summary["total"] == 0
        assert summary["by_category"] == {}

    def test_missing_project(self, reg):
        summary = reg.project_data_summary("nonexistent")
        assert summary["total"] == 0

    def test_with_mixed_data(self, reg):
        reg.add_project("proj", platform="visium")
        reg.register_biodata("proj", "spatial_seq", "visium", "/a.h5ad")
        reg.register_biodata("proj", "spatial_seq", "visium", "/b.h5ad")
        reg.register_biodata("proj", "image", "imc", "/c.tiff", fmt="tiff")
        summary = reg.project_data_summary("proj")
        assert summary["total"] == 3
        assert summary["by_category"]["spatial_seq"] == 2
        assert summary["by_category"]["image"] == 1
        assert summary["by_platform"]["visium"] == 2
        assert summary["by_platform"]["imc"] == 1


class TestMigrateDatasetsToBioData:
    """Migration tests. Note: register_dataset now dual-writes to BioData,
    so these tests verify that migrate_datasets_to_biodata() correctly
    handles already-migrated (dual-written) datasets and legacy datasets
    without bio_data_id."""

    def test_migrate_already_dual_written_is_noop(self, reg):
        """Datasets created by dual-write already have bio_data_id, so migration skips them."""
        reg.add_project("proj", platform="visium")
        reg.register_dataset("proj", phase="qc_filter", uri="/a.h5ad")
        reg.register_dataset("proj", phase="preprocess", uri="/b.h5ad")
        count = reg.migrate_datasets_to_biodata()
        assert count == 0  # dual-write already handled these
        biodata = reg.list_biodata(project_name="proj")
        assert len(biodata) == 2
        assert all(bd["type"] == "spatial_seq" for bd in biodata)

    def test_migrate_legacy_dataset_without_biodata(self, reg):
        """Simulate a legacy dataset (pre-dual-write) by clearing bio_data_id."""
        reg.add_project("proj", platform="imc")
        ds_id = reg.register_dataset(
            "proj", phase="ingest_load", uri="/img.tiff", fmt="tiff", file_role="image"
        )
        # Remove the dual-written BioData to simulate legacy data
        with reg._session() as sess:
            ds = sess.get(reg._Dataset, ds_id)
            if ds.bio_data_id is not None:
                bd = sess.get(reg._BioData, ds.bio_data_id)
                if bd is not None:
                    sess.delete(bd)
                ds.bio_data_id = None
                sess.commit()
        count = reg.migrate_datasets_to_biodata()
        assert count == 1

    def test_migrate_preserves_legacy_link(self, reg):
        """Legacy migration sets legacy_dataset_id on the new BioData row."""
        reg.add_project("proj", platform="visium")
        reg.register_dataset("proj", phase="qc_filter", uri="/a.h5ad")
        # Dual-write already created a BioData — verify linkage
        biodata = reg.list_biodata(project_name="proj")
        assert len(biodata) == 1
        # Dual-write does NOT set legacy_dataset_id (only migration does)
        # But the forward link on dataset IS set
        datasets = reg.list_datasets(project_name="proj")
        assert datasets[0].get("bio_data_id") is not None

    def test_migrate_idempotent(self, reg):
        reg.add_project("proj", platform="visium")
        reg.register_dataset("proj", phase="qc_filter", uri="/a.h5ad")
        # Dual-write already handled it, so both calls return 0
        count1 = reg.migrate_datasets_to_biodata()
        count2 = reg.migrate_datasets_to_biodata()
        assert count1 == 0
        assert count2 == 0

    def test_dual_write_imc_h5ad_uses_platform_registry(self, reg):
        """IMC h5ad files should be classified via platform registry as BioImage."""
        reg.add_project("proj", platform="imc")
        reg.register_dataset("proj", phase="qc_filter", uri="/results/adata.filtered.h5ad")
        biodata = reg.list_biodata(project_name="proj")
        assert len(biodata) == 1
        assert biodata[0]["type"] == "image"

    def test_dual_write_chromium_uses_platform_registry(self, reg):
        """chromium_3p datasets should be classified as rnaseq via platform registry."""
        reg.add_project("proj", platform="chromium_3p")
        reg.register_dataset("proj", phase="qc_filter", uri="/results/counts.h5ad")
        biodata = reg.list_biodata(project_name="proj")
        assert len(biodata) == 1
        assert biodata[0]["type"] == "rnaseq"


class TestDeprecationAndWarnings:
    def test_register_dataset_emits_deprecation(self, reg):
        reg.add_project("proj", platform="visium")
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            reg.register_dataset("proj", phase="qc_filter", uri="/a.h5ad")
            deprecation_warns = [x for x in w if issubclass(x.category, DeprecationWarning)]
            assert len(deprecation_warns) == 1
            assert "register_biodata" in str(deprecation_warns[0].message)

    def test_register_biodata_ignores_unknown_kwargs(self, reg):
        """Unknown type_kwargs are silently dropped (with a log warning)."""
        reg.add_project("proj", platform="visium")
        bd_id = reg.register_biodata(
            "proj",
            "spatial_seq",
            "visium",
            "/a.h5ad",
            totally_bogus_field="xyz",
        )
        assert isinstance(bd_id, int)
        # Verify the data was still created successfully
        bd = reg.get_biodata(bd_id)
        assert bd is not None
        assert bd["platform"] == "visium"


class TestStatusIncludesNewCounts:
    def test_status_has_subject_sample_biodata_counts(self, reg):
        s = reg.status()
        assert "n_subjects" in s
        assert "n_samples" in s
        assert "n_biodata" in s
        assert s["n_subjects"] == 0
        assert s["n_samples"] == 0
        assert s["n_biodata"] == 0

    def test_status_counts_after_adds(self, reg):
        reg.add_project("proj", platform="visium")
        reg.add_subject("PT001")
        reg.add_sample("S001", subject_id="PT001", project_name="proj")
        reg.register_biodata("proj", "spatial_seq", "visium", "/a.h5ad")
        s = reg.status()
        assert s["n_subjects"] == 1
        assert s["n_samples"] == 1
        assert s["n_biodata"] == 1


# ---------------------------------------------------------------------------
# Modality auto-fill
# ---------------------------------------------------------------------------


class TestModalityAutoFill:
    def test_modality_auto_filled_from_platform(self, reg):
        reg.add_project("proj", platform="imc")
        bd_id = reg.register_biodata("proj", "image", "imc", "/data/x.tiff", fmt="tiff")
        bd = reg.get_biodata(bd_id)
        assert bd["modality"] == "Spatial Proteomics - Mass Spec"

    def test_modality_auto_filled_visium(self, reg):
        reg.add_project("proj", platform="visium")
        bd_id = reg.register_biodata("proj", "spatial_seq", "visium", "/a.h5ad")
        bd = reg.get_biodata(bd_id)
        assert bd["modality"] == "Spatial Transcriptomics - Sequencing"

    def test_modality_auto_filled_xenium(self, reg):
        reg.add_project("proj", platform="xenium")
        bd_id = reg.register_biodata("proj", "spatial_seq", "xenium", "/a.h5ad")
        bd = reg.get_biodata(bd_id)
        assert bd["modality"] == "Spatial Transcriptomics - ISH"

    def test_modality_none_for_unknown_platform(self, reg):
        reg.add_project("proj", platform="visium")
        bd_id = reg.register_biodata("proj", "spatial_seq", "unknown_plat", "/a.h5ad")
        bd = reg.get_biodata(bd_id)
        assert bd["modality"] is None


# ---------------------------------------------------------------------------
# Defaults auto-fill from PlatformSpec
# ---------------------------------------------------------------------------


class TestDefaultsAutoFill:
    def test_imc_autofill_staining_protocol(self, reg):
        """IMC should auto-fill staining_protocol and image_type from defaults."""
        reg.add_project("proj", platform="imc")
        bd_id = reg.register_biodata("proj", "image", "imc", "/data/x.tiff", fmt="tiff")
        bd = reg.get_biodata(bd_id)
        assert bd["staining_protocol"] == "IMC"
        assert bd["image_type"] == "multiplexed"

    def test_explicit_overrides_default(self, reg):
        """Explicit kwargs should override PlatformSpec defaults."""
        reg.add_project("proj", platform="imc")
        bd_id = reg.register_biodata(
            "proj", "image", "imc", "/data/x.tiff", fmt="tiff", staining_protocol="custom_protocol"
        )
        bd = reg.get_biodata(bd_id)
        assert bd["staining_protocol"] == "custom_protocol"

    def test_visium_hd_autofill_bin_size(self, reg):
        reg.add_project("proj", platform="visium_hd")
        bd_id = reg.register_biodata("proj", "spatial_seq", "visium_hd", "/a.h5ad")
        bd = reg.get_biodata(bd_id)
        assert bd["bin_size_um"] == 8.0
        assert bd["spatial_resolution"] == "spot"

    def test_xenium_autofill_coordinate_system(self, reg):
        reg.add_project("proj", platform="xenium")
        bd_id = reg.register_biodata("proj", "spatial_seq", "xenium", "/a.h5ad")
        bd = reg.get_biodata(bd_id)
        assert bd["spatial_resolution"] == "single_cell"
        assert bd["coordinate_system"] == "micron"

    def test_chromium_autofill_chemistry(self, reg):
        reg.add_project("proj", platform="chromium_3p")
        bd_id = reg.register_biodata("proj", "rnaseq", "chromium_3p", "/a.h5ad")
        bd = reg.get_biodata(bd_id)
        assert bd["chemistry"] == "chromium_v3"
        assert bd["library_type"] == "single_cell"


# ---------------------------------------------------------------------------
# Format inference from URI
# ---------------------------------------------------------------------------


class TestFormatInference:
    def test_h5ad_inferred(self, reg):
        reg.add_project("proj", platform="visium")
        bd_id = reg.register_biodata("proj", "spatial_seq", "visium", "/data/adata.h5ad")
        bd = reg.get_biodata(bd_id)
        assert bd["format"] == "h5ad"

    def test_tiff_inferred(self, reg):
        reg.add_project("proj", platform="imc")
        bd_id = reg.register_biodata("proj", "image", "imc", "/data/image.tiff")
        bd = reg.get_biodata(bd_id)
        assert bd["format"] == "tiff"

    def test_zarr_inferred(self, reg):
        reg.add_project("proj", platform="visium_hd")
        bd_id = reg.register_biodata("proj", "spatial_seq", "visium_hd", "/data/sdata.zarr")
        bd = reg.get_biodata(bd_id)
        assert bd["format"] == "zarr"

    def test_explicit_format_not_overridden(self, reg):
        reg.add_project("proj", platform="visium")
        bd_id = reg.register_biodata(
            "proj", "spatial_seq", "visium", "/data/adata.h5ad", fmt="custom"
        )
        bd = reg.get_biodata(bd_id)
        assert bd["format"] == "custom"

    def test_no_extension_stays_none(self, reg):
        reg.add_project("proj", platform="visium")
        bd_id = reg.register_biodata("proj", "spatial_seq", "visium", "/data/noext")
        bd = reg.get_biodata(bd_id)
        assert bd["format"] is None


# ---------------------------------------------------------------------------
# Dual-write: register_dataset -> BioData
# ---------------------------------------------------------------------------


class TestDualWrite:
    def test_register_dataset_creates_biodata(self, reg):
        reg.add_project("proj", platform="visium")
        reg.register_dataset("proj", phase="qc_filter", uri="/a.h5ad")
        # Should have created a BioData row
        biodata = reg.list_biodata(project_name="proj")
        assert len(biodata) == 1
        assert biodata[0]["platform"] == "visium"
        assert biodata[0]["type"] == "spatial_seq"

    def test_dual_write_links_forward(self, reg):
        reg.add_project("proj", platform="imc")
        reg.register_dataset("proj", phase="qc_filter", uri="/a.h5ad")
        # Check that dataset has bio_data_id
        datasets = reg.list_datasets(project_name="proj")
        assert len(datasets) == 1
        assert datasets[0].get("bio_data_id") is not None

    def test_dual_write_preserves_metadata(self, reg):
        reg.add_project("proj", platform="visium")
        reg.register_dataset(
            "proj",
            phase="scoring",
            uri="/scored.h5ad",
            n_obs=50000,
            n_vars=3000,
            file_role="primary",
        )
        biodata = reg.list_biodata(project_name="proj")
        assert len(biodata) == 1
        assert biodata[0]["n_obs"] == 50000
        assert biodata[0]["n_vars"] == 3000
        assert biodata[0]["phase"] == "scoring"

    def test_dual_write_imc_creates_image_type(self, reg):
        reg.add_project("proj", platform="imc")
        reg.register_dataset("proj", phase="ingest_load", uri="/img.h5ad")
        biodata = reg.list_biodata(project_name="proj")
        assert len(biodata) == 1
        assert biodata[0]["type"] == "image"
        assert biodata[0]["modality"] == "Spatial Proteomics - Mass Spec"
