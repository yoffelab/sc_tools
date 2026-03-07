"""Unit tests for sc_tools.registry."""

from __future__ import annotations

import json

import pytest

# Skip all tests if SQLAlchemy is not installed
sqlalchemy = pytest.importorskip("sqlalchemy", reason="sqlalchemy not installed")


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
