"""Unit tests for sc_tools.registry -- four-layer schema."""

from __future__ import annotations

import warnings

import pytest

# Skip all tests if SQLAlchemy is not installed
sqlalchemy = pytest.importorskip("sqlalchemy", reason="sqlalchemy not installed")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def reg(tmp_path):
    """Return a file-backed SQLite Registry for testing."""
    from sc_tools.registry import Registry

    db_url = f"sqlite:///{tmp_path / 'test_registry.db'}"
    return Registry(db_url=db_url)


def _make_data_source(reg, name="src_geo_ibd", **kwargs):
    """Helper: register a data source with sensible defaults."""
    defaults = {
        "uri": "s3://bucket/raw/geo_ibd",
        "platform": "visium",
        "source_type": "geo",
        "disease": "ibd",
        "tissue": "colon",
    }
    defaults.update(kwargs)
    return reg.register_data_source(name, **defaults)


def _make_inventory_item(reg, name="item_rna_v1", **kwargs):
    """Helper: register an inventory item with sensible defaults."""
    defaults = {"uri": "/data/rna_v1.h5ad", "modality": "rna"}
    defaults.update(kwargs)
    return reg.register_inventory_item(name, **defaults)


def _make_dataset(reg, name="ds_ibd_v1"):
    """Helper: create a dataset."""
    return reg.create_dataset(name)


def _setup_project_with_dataset(reg, project="proj_a", dataset="ds_a"):
    """Helper: create a project, dataset, and link them."""
    reg.add_project(project, platform="visium")
    reg.create_dataset(dataset)
    reg.link_project_dataset(project, dataset)
    return project, dataset


# ---------------------------------------------------------------------------
# TestProjects
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
        reg.add_project("proj_b")
        proj = reg.get_project("proj_b")
        assert proj is not None
        assert proj["name"] == "proj_b"

    def test_get_project_not_found(self, reg):
        assert reg.get_project("nonexistent") is None

    def test_list_projects(self, reg):
        reg.add_project("proj_a")
        reg.add_project("proj_b")
        projects = reg.list_projects()
        names = [p["name"] for p in projects]
        assert "proj_a" in names
        assert "proj_b" in names

    def test_delete_project(self, reg):
        reg.add_project("proj_del")
        assert reg.delete_project("proj_del") is True
        assert reg.get_project("proj_del") is None

    def test_delete_project_not_found(self, reg):
        assert reg.delete_project("nonexistent") is False


# ---------------------------------------------------------------------------
# TestDataSources
# ---------------------------------------------------------------------------


class TestDataSources:
    def test_register_data_source(self, reg):
        sid = _make_data_source(reg)
        assert isinstance(sid, int)
        assert sid > 0

    def test_register_data_source_idempotent(self, reg):
        sid1 = _make_data_source(reg)
        sid2 = _make_data_source(reg)
        assert sid1 == sid2

    def test_get_data_source(self, reg):
        _make_data_source(reg, name="src_a")
        ds = reg.get_data_source("src_a")
        assert ds is not None
        assert ds["name"] == "src_a"
        assert ds["platform"] == "visium"

    def test_get_data_source_not_found(self, reg):
        assert reg.get_data_source("nonexistent") is None

    def test_list_data_sources_no_filter(self, reg):
        _make_data_source(reg, name="src_a", platform="visium")
        _make_data_source(reg, name="src_b", platform="imc")
        results = reg.list_data_sources()
        assert len(results) == 2

    def test_list_data_sources_filter_platform(self, reg):
        _make_data_source(reg, name="src_a", platform="visium")
        _make_data_source(reg, name="src_b", platform="imc")
        results = reg.list_data_sources(platform="visium")
        assert len(results) == 1
        assert results[0]["name"] == "src_a"

    def test_list_data_sources_filter_disease(self, reg):
        _make_data_source(reg, name="src_a", disease="ibd")
        _make_data_source(reg, name="src_b", disease="cancer")
        results = reg.list_data_sources(disease="ibd")
        assert len(results) == 1
        assert results[0]["name"] == "src_a"

    def test_list_data_sources_filter_tissue(self, reg):
        _make_data_source(reg, name="src_a", tissue="colon")
        _make_data_source(reg, name="src_b", tissue="liver")
        results = reg.list_data_sources(tissue="colon")
        assert len(results) == 1
        assert results[0]["name"] == "src_a"

    def test_list_data_sources_filter_source_type(self, reg):
        _make_data_source(reg, name="src_a", source_type="geo")
        _make_data_source(reg, name="src_b", source_type="hpc_directory")
        results = reg.list_data_sources(source_type="geo")
        assert len(results) == 1
        assert results[0]["name"] == "src_a"


# ---------------------------------------------------------------------------
# TestInventoryItems
# ---------------------------------------------------------------------------


class TestInventoryItems:
    def test_register_linked_to_data_source(self, reg):
        _make_data_source(reg, name="src_a")
        iid = reg.register_inventory_item(
            "item_a",
            uri="/data/item_a.h5ad",
            modality="rna",
            data_source_name="src_a",
        )
        assert isinstance(iid, int)
        item = reg.get_inventory_item("item_a")
        assert item is not None
        assert item["data_source_id"] is not None

    def test_register_standalone(self, reg):
        iid = _make_inventory_item(reg, name="item_standalone")
        assert isinstance(iid, int)
        item = reg.get_inventory_item("item_standalone")
        assert item is not None
        assert item["data_source_id"] is None

    def test_register_idempotent(self, reg):
        iid1 = _make_inventory_item(reg, name="item_x")
        iid2 = _make_inventory_item(reg, name="item_x")
        assert iid1 == iid2

    def test_register_invalid_data_source_raises(self, reg):
        with pytest.raises(ValueError, match="DataSource.*not found"):
            reg.register_inventory_item(
                "item_bad",
                uri="/data/bad.h5ad",
                modality="rna",
                data_source_name="nonexistent_source",
            )

    def test_get_inventory_item(self, reg):
        _make_inventory_item(reg, name="item_a")
        item = reg.get_inventory_item("item_a")
        assert item["name"] == "item_a"
        assert item["modality"] == "rna"

    def test_get_inventory_item_not_found(self, reg):
        assert reg.get_inventory_item("nonexistent") is None

    def test_list_inventory_items_no_filter(self, reg):
        _make_inventory_item(reg, name="item_a", modality="rna")
        _make_inventory_item(reg, name="item_b", modality="protein")
        results = reg.list_inventory_items()
        assert len(results) == 2

    def test_list_inventory_items_filter_modality(self, reg):
        _make_inventory_item(reg, name="item_a", modality="rna")
        _make_inventory_item(reg, name="item_b", modality="protein")
        results = reg.list_inventory_items(modality="rna")
        assert len(results) == 1
        assert results[0]["name"] == "item_a"

    def test_list_inventory_items_filter_platform(self, reg):
        _make_inventory_item(reg, name="item_a", platform="visium")
        _make_inventory_item(reg, name="item_b", platform="imc")
        results = reg.list_inventory_items(platform="visium")
        assert len(results) == 1
        assert results[0]["name"] == "item_a"


# ---------------------------------------------------------------------------
# TestDatasets
# ---------------------------------------------------------------------------


class TestDatasets:
    def test_create_dataset(self, reg):
        did = reg.create_dataset("ds_a")
        assert isinstance(did, int)
        assert did > 0

    def test_create_dataset_idempotent(self, reg):
        did1 = reg.create_dataset("ds_a")
        did2 = reg.create_dataset("ds_a")
        assert did1 == did2

    def test_get_dataset_by_name(self, reg):
        reg.create_dataset("ds_a", description="test dataset")
        ds = reg.get_dataset("ds_a")
        assert ds is not None
        assert ds["name"] == "ds_a"
        assert ds["version"] == 1
        assert ds["is_current"] is True
        assert ds["description"] == "test dataset"

    def test_get_dataset_by_version(self, reg):
        reg.create_dataset("ds_a")
        ds = reg.get_dataset("ds_a", version=1)
        assert ds is not None
        assert ds["version"] == 1

    def test_get_dataset_not_found(self, reg):
        assert reg.get_dataset("nonexistent") is None

    def test_get_dataset_current(self, reg):
        reg.create_dataset("ds_a")
        ds = reg.get_dataset("ds_a")
        assert ds["is_current"] is True

    def test_add_dataset_member(self, reg):
        reg.create_dataset("ds_a")
        _make_inventory_item(reg, name="item_rna")
        mid = reg.add_dataset_member("ds_a", "item_rna", "rna")
        assert isinstance(mid, int)

    def test_add_dataset_member_nonexistent_dataset_raises(self, reg):
        _make_inventory_item(reg, name="item_rna")
        with pytest.raises(ValueError, match="Dataset.*not found"):
            reg.add_dataset_member("nonexistent_ds", "item_rna", "rna")

    def test_add_dataset_member_nonexistent_inventory_raises(self, reg):
        reg.create_dataset("ds_a")
        with pytest.raises(ValueError, match="Inventory item.*not found"):
            reg.add_dataset_member("ds_a", "nonexistent_item", "rna")

    def test_add_dataset_member_duplicate_modality_key_raises(self, reg):
        reg.create_dataset("ds_a")
        _make_inventory_item(reg, name="item_rna")
        _make_inventory_item(reg, name="item_rna2", uri="/data/rna2.h5ad")
        reg.add_dataset_member("ds_a", "item_rna", "rna")
        with pytest.raises(ValueError, match="Modality key.*already exists"):
            reg.add_dataset_member("ds_a", "item_rna2", "rna")

    def test_remove_dataset_member(self, reg):
        reg.create_dataset("ds_a")
        _make_inventory_item(reg, name="item_rna")
        reg.add_dataset_member("ds_a", "item_rna", "rna")
        reg.remove_dataset_member("ds_a", "rna")
        members = reg.get_dataset_members("ds_a")
        assert len(members) == 0

    def test_remove_dataset_member_nonexistent_dataset_raises(self, reg):
        with pytest.raises(ValueError, match="Dataset.*not found"):
            reg.remove_dataset_member("nonexistent_ds", "rna")

    def test_remove_dataset_member_nonexistent_key_raises(self, reg):
        reg.create_dataset("ds_a")
        with pytest.raises(ValueError, match="Modality key.*not found"):
            reg.remove_dataset_member("ds_a", "nonexistent_key")

    def test_get_dataset_members_enriched(self, reg):
        reg.create_dataset("ds_a")
        _make_inventory_item(reg, name="item_rna", platform="visium", n_obs=5000, n_vars=2000)
        reg.add_dataset_member("ds_a", "item_rna", "rna")
        members = reg.get_dataset_members("ds_a")
        assert len(members) == 1
        m = members[0]
        assert m["modality_key"] == "rna"
        assert m["inventory_name"] == "item_rna"
        assert m["inventory_uri"] == "/data/rna_v1.h5ad"
        assert m["inventory_platform"] == "visium"
        assert m["n_obs"] == 5000
        assert m["n_vars"] == 2000

    def test_list_datasets_all(self, reg):
        reg.create_dataset("ds_a")
        reg.create_dataset("ds_b")
        results = reg.list_datasets()
        names = [d["name"] for d in results]
        assert "ds_a" in names
        assert "ds_b" in names

    def test_list_datasets_by_project(self, reg):
        _setup_project_with_dataset(reg, "proj_a", "ds_a")
        reg.create_dataset("ds_b")  # not linked to any project
        results = reg.list_datasets(project_name="proj_a")
        assert len(results) == 1
        assert results[0]["name"] == "ds_a"


# ---------------------------------------------------------------------------
# TestDatasetVersioning
# ---------------------------------------------------------------------------


class TestDatasetVersioning:
    def test_bump_version_copies_members(self, reg):
        reg.create_dataset("ds_a")
        _make_inventory_item(reg, name="item_rna")
        _make_inventory_item(reg, name="item_protein", uri="/data/protein.h5ad", modality="protein")
        reg.add_dataset_member("ds_a", "item_rna", "rna")
        reg.add_dataset_member("ds_a", "item_protein", "protein")

        new_id = reg.bump_dataset_version("ds_a")
        assert isinstance(new_id, int)

        members_v2 = reg.get_dataset_members("ds_a")
        keys = {m["modality_key"] for m in members_v2}
        assert keys == {"rna", "protein"}

    def test_old_version_not_current(self, reg):
        reg.create_dataset("ds_a")
        reg.bump_dataset_version("ds_a")
        old = reg.get_dataset("ds_a", version=1)
        assert old["is_current"] is False

    def test_new_version_is_current(self, reg):
        reg.create_dataset("ds_a")
        reg.bump_dataset_version("ds_a")
        new = reg.get_dataset("ds_a", version=2)
        assert new["is_current"] is True

    def test_get_dataset_returns_current(self, reg):
        reg.create_dataset("ds_a")
        reg.bump_dataset_version("ds_a")
        ds = reg.get_dataset("ds_a")
        assert ds["version"] == 2
        assert ds["is_current"] is True

    def test_get_dataset_old_version(self, reg):
        reg.create_dataset("ds_a")
        reg.bump_dataset_version("ds_a")
        ds = reg.get_dataset("ds_a", version=1)
        assert ds["version"] == 1

    def test_bump_nonexistent_raises(self, reg):
        with pytest.raises(ValueError, match="Dataset.*not found"):
            reg.bump_dataset_version("nonexistent")


# ---------------------------------------------------------------------------
# TestProjectDatasets
# ---------------------------------------------------------------------------


class TestProjectDatasets:
    def test_link_project_dataset(self, reg):
        reg.add_project("proj_a", platform="visium")
        reg.create_dataset("ds_a")
        lid = reg.link_project_dataset("proj_a", "ds_a")
        assert isinstance(lid, int)

    def test_link_project_dataset_idempotent(self, reg):
        reg.add_project("proj_a", platform="visium")
        reg.create_dataset("ds_a")
        lid1 = reg.link_project_dataset("proj_a", "ds_a")
        lid2 = reg.link_project_dataset("proj_a", "ds_a")
        assert lid1 == lid2

    def test_link_project_dataset_nonexistent_project_raises(self, reg):
        reg.create_dataset("ds_a")
        with pytest.raises(ValueError, match="Project.*not found"):
            reg.link_project_dataset("nonexistent", "ds_a")

    def test_link_project_dataset_nonexistent_dataset_raises(self, reg):
        reg.add_project("proj_a")
        with pytest.raises(ValueError, match="Dataset.*not found"):
            reg.link_project_dataset("proj_a", "nonexistent")

    def test_list_project_datasets(self, reg):
        reg.add_project("proj_a", platform="visium")
        reg.create_dataset("ds_a")
        reg.create_dataset("ds_b")
        reg.link_project_dataset("proj_a", "ds_a")
        reg.link_project_dataset("proj_a", "ds_b", role="validation")
        results = reg.list_project_datasets("proj_a")
        assert len(results) == 2
        roles = {d["link_role"] for d in results}
        assert "primary" in roles
        assert "validation" in roles

    def test_link_different_roles(self, reg):
        reg.add_project("proj_a")
        reg.create_dataset("ds_primary")
        reg.create_dataset("ds_validation")
        reg.create_dataset("ds_reference")
        reg.link_project_dataset("proj_a", "ds_primary", role="primary")
        reg.link_project_dataset("proj_a", "ds_validation", role="validation")
        reg.link_project_dataset("proj_a", "ds_reference", role="reference")
        results = reg.list_project_datasets("proj_a")
        assert len(results) == 3


# ---------------------------------------------------------------------------
# TestProjectPhases
# ---------------------------------------------------------------------------


class TestProjectPhases:
    def test_upsert_phase_data_processing(self, reg):
        proj, ds = _setup_project_with_dataset(reg)
        reg.upsert_phase(proj, ds, "data_processing", "qc_filter", status="in_progress")
        phase = reg.get_phase(proj, "data_processing", "qc_filter")
        assert phase is not None
        assert phase["status"] == "in_progress"
        assert phase["phase_group"] == "data_processing"
        assert phase["subphase"] == "qc_filter"

    def test_upsert_phase_discovery(self, reg):
        proj, ds = _setup_project_with_dataset(reg)
        reg.upsert_phase(proj, ds, "discovery", "clustering_v1", status="in_progress")
        phase = reg.get_phase(proj, "discovery", "clustering_v1")
        assert phase is not None
        assert phase["phase_group"] == "discovery"
        assert phase["subphase"] == "clustering_v1"

    def test_upsert_phase_update(self, reg):
        proj, ds = _setup_project_with_dataset(reg)
        reg.upsert_phase(proj, ds, "data_processing", "qc_filter", status="in_progress")
        reg.upsert_phase(
            proj,
            ds,
            "data_processing",
            "qc_filter",
            status="ready",
            n_obs=5000,
            n_vars=2000,
        )
        phase = reg.get_phase(proj, "data_processing", "qc_filter")
        assert phase["status"] == "ready"
        assert phase["n_obs"] == 5000
        assert phase["n_vars"] == 2000

    def test_get_phase(self, reg):
        proj, ds = _setup_project_with_dataset(reg)
        reg.upsert_phase(proj, ds, "data_processing", "preprocess", status="in_progress")
        phase = reg.get_phase(proj, "data_processing", "preprocess")
        assert phase is not None
        assert phase["subphase"] == "preprocess"

    def test_get_phase_not_found(self, reg):
        reg.add_project("proj_a")
        assert reg.get_phase("proj_a", "data_processing", "nonexistent") is None

    def test_list_phases_all(self, reg):
        proj, ds = _setup_project_with_dataset(reg)
        reg.upsert_phase(proj, ds, "data_processing", "qc_filter")
        reg.upsert_phase(proj, ds, "data_processing", "preprocess")
        reg.upsert_phase(proj, ds, "discovery", "clustering_v1")
        phases = reg.list_phases(proj)
        assert len(phases) == 3

    def test_list_phases_filtered_by_group(self, reg):
        proj, ds = _setup_project_with_dataset(reg)
        reg.upsert_phase(proj, ds, "data_processing", "qc_filter")
        reg.upsert_phase(proj, ds, "data_processing", "preprocess")
        reg.upsert_phase(proj, ds, "discovery", "clustering_v1")
        dp_phases = reg.list_phases(proj, phase_group="data_processing")
        assert len(dp_phases) == 2
        disc_phases = reg.list_phases(proj, phase_group="discovery")
        assert len(disc_phases) == 1

    def test_mark_phase_complete(self, reg):
        proj, ds = _setup_project_with_dataset(reg)
        reg.upsert_phase(proj, ds, "data_processing", "qc_filter", status="in_progress")
        reg.mark_phase_complete(proj, "data_processing", "qc_filter")
        phase = reg.get_phase(proj, "data_processing", "qc_filter")
        assert phase["status"] == "ready"

    def test_mark_phase_complete_not_found_raises(self, reg):
        reg.add_project("proj_a")
        with pytest.raises(ValueError, match="Phase.*not found"):
            reg.mark_phase_complete("proj_a", "data_processing", "nonexistent")

    def test_upsert_phase_nonexistent_project_raises(self, reg):
        reg.create_dataset("ds_a")
        with pytest.raises(ValueError, match="Project.*not found"):
            reg.upsert_phase("nonexistent", "ds_a", "data_processing", "qc_filter")

    def test_upsert_phase_nonexistent_dataset_raises(self, reg):
        reg.add_project("proj_a")
        with pytest.raises(ValueError, match="Dataset.*not found"):
            reg.upsert_phase("proj_a", "nonexistent", "data_processing", "qc_filter")


# ---------------------------------------------------------------------------
# TestPatients
# ---------------------------------------------------------------------------


class TestPatients:
    def test_add_patient(self, reg):
        pid = reg.add_patient("PT001", metadata={"organism": "human", "sex": "M"})
        assert isinstance(pid, int)
        assert pid > 0

    def test_add_patient_idempotent(self, reg):
        pid1 = reg.add_patient("PT001")
        pid2 = reg.add_patient("PT001")
        assert pid1 == pid2

    def test_get_patient(self, reg):
        reg.add_patient("PT001", metadata={"organism": "human"})
        pat = reg.get_patient("PT001")
        assert pat is not None
        assert pat["patient_id"] == "PT001"
        assert pat["meta"]["organism"] == "human"

    def test_get_patient_not_found(self, reg):
        assert reg.get_patient("NONEXISTENT") is None

    def test_list_patients(self, reg):
        reg.add_patient("PT001")
        reg.add_patient("PT002")
        patients = reg.list_patients()
        ids = [p["patient_id"] for p in patients]
        assert "PT001" in ids
        assert "PT002" in ids


# ---------------------------------------------------------------------------
# TestSamples
# ---------------------------------------------------------------------------


class TestSamples:
    def test_register_sample(self, reg):
        sid = reg.register_sample("SAMPLE001", tissue="colon")
        assert isinstance(sid, int)
        assert sid > 0

    def test_register_sample_idempotent(self, reg):
        sid1 = reg.register_sample("SAMPLE001")
        sid2 = reg.register_sample("SAMPLE001")
        assert sid1 == sid2

    def test_register_sample_with_patient(self, reg):
        reg.add_patient("PT001")
        reg.register_sample("PT001_biopsy1", patient_id="PT001", tissue="colon")
        sample = reg.get_sample("PT001_biopsy1")
        assert sample is not None
        assert sample["sample_id"] == "PT001_biopsy1"
        assert sample["tissue"] == "colon"
        assert sample["patient_id"] is not None

    def test_register_sample_patient_not_found(self, reg):
        with pytest.raises(ValueError, match="Patient.*not found"):
            reg.register_sample("SAMPLE001", patient_id="NONEXISTENT")

    def test_get_sample(self, reg):
        reg.register_sample("SAMPLE001", tissue="liver")
        s = reg.get_sample("SAMPLE001")
        assert s is not None
        assert s["tissue"] == "liver"

    def test_get_sample_not_found(self, reg):
        assert reg.get_sample("NONEXISTENT") is None

    def test_list_samples(self, reg):
        reg.register_sample("S001")
        reg.register_sample("S002")
        samples = reg.list_samples()
        ids = [s["sample_id"] for s in samples]
        assert "S001" in ids
        assert "S002" in ids

    def test_list_samples_by_patient(self, reg):
        reg.add_patient("PT001")
        reg.add_patient("PT002")
        reg.register_sample("S001", patient_id="PT001")
        reg.register_sample("S002", patient_id="PT002")
        samples = reg.list_samples(patient_id="PT001")
        assert len(samples) == 1
        assert samples[0]["sample_id"] == "S001"

    def test_list_samples_by_tissue(self, reg):
        reg.register_sample("S001", tissue="colon")
        reg.register_sample("S002", tissue="liver")
        samples = reg.list_samples(tissue="colon")
        assert len(samples) == 1
        assert samples[0]["sample_id"] == "S001"

    def test_inventory_item_with_sample(self, reg):
        reg.add_patient("PT001")
        reg.register_sample("PT001_bx1", patient_id="PT001")
        reg.register_inventory_item(
            "rna_pt001",
            uri="/data/rna.h5ad",
            modality="rna",
            sample_name="PT001_bx1",
        )
        item = reg.get_inventory_item("rna_pt001")
        assert item["sample_id"] is not None

    def test_inventory_item_sample_not_found(self, reg):
        with pytest.raises(ValueError, match="Sample.*not found"):
            reg.register_inventory_item(
                "rna_bad",
                uri="/data/rna.h5ad",
                modality="rna",
                sample_name="NONEXISTENT",
            )

    def test_patient_to_dataset_via_sample(self, reg):
        """Full chain: patient -> sample -> inventory -> dataset -> project."""
        reg.add_project("proj_a")
        reg.add_patient("PT001", metadata={"organism": "human"})
        reg.register_sample("PT001_bx1", patient_id="PT001", tissue="colon")
        reg.register_inventory_item(
            "rna_pt001",
            uri="/data/rna.h5ad",
            modality="rna",
            sample_name="PT001_bx1",
        )
        reg.create_dataset("ds_a")
        reg.add_dataset_member("ds_a", "rna_pt001", "rna")
        reg.link_project_dataset("proj_a", "ds_a")

        # list_subjects with project filter should find PT001
        subjects = reg.list_subjects(project_name="proj_a")
        assert len(subjects) == 1
        assert subjects[0]["patient_id"] == "PT001"


# ---------------------------------------------------------------------------
# TestProvenance
# ---------------------------------------------------------------------------


class TestProvenance:
    def test_record_provenance_for_phase(self, reg):
        proj, ds = _setup_project_with_dataset(reg)
        reg.upsert_phase(proj, ds, "data_processing", "qc_filter")
        with reg._session() as sess:
            phase_row = (
                sess.query(reg._ProjectPhase)
                .filter_by(phase_group="data_processing", subphase="qc_filter")
                .first()
            )
            phase_id = phase_row.id

        prov_id = reg.record_provenance(
            "sc_tools.qc",
            tool_version="1.10.0",
            phase_id=phase_id,
            n_input_obs=10000,
            n_output_obs=8000,
        )
        assert isinstance(prov_id, int)

    def test_record_provenance_for_dataset(self, reg):
        ds_id = reg.create_dataset("ds_prov")
        prov_id = reg.record_provenance("mudata_builder", dataset_id=ds_id)
        assert isinstance(prov_id, int)

    def test_get_provenance_by_phase_id(self, reg):
        proj, ds = _setup_project_with_dataset(reg)
        reg.upsert_phase(proj, ds, "data_processing", "qc_filter")
        with reg._session() as sess:
            phase_row = (
                sess.query(reg._ProjectPhase)
                .filter_by(phase_group="data_processing", subphase="qc_filter")
                .first()
            )
            phase_id = phase_row.id
        reg.record_provenance("sc_tools.qc", phase_id=phase_id)
        reg.record_provenance("scvi-tools", phase_id=phase_id)
        records = reg.get_provenance(phase_id=phase_id)
        assert len(records) == 2
        tools = {r["tool"] for r in records}
        assert tools == {"sc_tools.qc", "scvi-tools"}

    def test_get_provenance_by_dataset_id(self, reg):
        ds_id = reg.create_dataset("ds_prov")
        reg.record_provenance("muon.concat", dataset_id=ds_id)
        records = reg.get_provenance(dataset_id=ds_id)
        assert len(records) == 1
        assert records[0]["tool"] == "muon.concat"

    def test_record_provenance_no_target_raises(self, reg):
        with pytest.raises(ValueError, match="Exactly one"):
            reg.record_provenance("scanpy")

    def test_record_provenance_multiple_targets_raises(self, reg):
        proj, ds = _setup_project_with_dataset(reg)
        reg.upsert_phase(proj, ds, "data_processing", "qc_filter")
        with reg._session() as sess:
            phase_row = (
                sess.query(reg._ProjectPhase)
                .filter_by(phase_group="data_processing", subphase="qc_filter")
                .first()
            )
            phase_id = phase_row.id
        ds_id = reg.create_dataset("ds_multi")
        with pytest.raises(ValueError, match="Exactly one"):
            reg.record_provenance("scanpy", phase_id=phase_id, dataset_id=ds_id)

    def test_record_provenance_params_stored(self, reg):
        ds_id = reg.create_dataset("ds_params")
        reg.record_provenance(
            "scanpy",
            dataset_id=ds_id,
            params={"min_genes": 200, "max_genes": 5000},
            environment={"python": "3.11"},
        )
        records = reg.get_provenance(dataset_id=ds_id)
        assert records[0]["params"]["min_genes"] == 200
        assert records[0]["environment"]["python"] == "3.11"


# ---------------------------------------------------------------------------
# TestBackwardCompat
# ---------------------------------------------------------------------------


class TestBackwardCompat:
    def test_register_data_emits_deprecation(self, reg):
        reg.add_project("proj_a")
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            reg.register_data("proj_a", "qc_filter", "/data/test.h5ad")
            assert len(w) == 1
            assert issubclass(w[0].category, DeprecationWarning)
            assert "register_data" in str(w[0].message)

    def test_register_dataset_emits_deprecation(self, reg):
        reg.add_project("proj_a")
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            reg.register_dataset("proj_a", "qc_filter", "/data/test.h5ad")
            assert len(w) == 1
            assert issubclass(w[0].category, DeprecationWarning)
            assert "register_dataset" in str(w[0].message)

    def test_register_biodata_emits_deprecation(self, reg):
        reg.add_project("proj_a")
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            reg.register_biodata("proj_a", "rna", "visium", "/data/test.h5ad")
            assert len(w) == 1
            assert issubclass(w[0].category, DeprecationWarning)
            assert "register_biodata" in str(w[0].message)

    def test_list_biodata_emits_deprecation(self, reg):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            reg.list_biodata()
            assert len(w) == 1
            assert issubclass(w[0].category, DeprecationWarning)
            assert "list_biodata" in str(w[0].message)

    def test_register_data_creates_inventory_item(self, reg):
        reg.add_project("proj_a")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            iid = reg.register_data("proj_a", "qc_filter", "/data/test.h5ad", category="rna")
        assert isinstance(iid, int)
        items = reg.list_inventory_items()
        assert len(items) == 1

    def test_register_biodata_creates_inventory_item(self, reg):
        reg.add_project("proj_a")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            iid = reg.register_biodata("proj_a", "rna", "visium", "/data/test.h5ad")
        assert isinstance(iid, int)
        items = reg.list_inventory_items()
        assert len(items) == 1
        assert items[0]["modality"] == "rna"


# ---------------------------------------------------------------------------
# TestStatus
# ---------------------------------------------------------------------------


class TestStatus:
    def test_status_counts(self, reg):
        reg.add_project("proj_a")
        _make_data_source(reg, name="src_a")
        _make_inventory_item(reg, name="item_a")
        reg.create_dataset("ds_a")
        reg.add_patient("PT001")
        reg.register_sample("S001", patient_id="PT001")

        s = reg.status()
        assert s["n_projects"] == 1
        assert s["n_data_sources"] == 1
        assert s["n_inventory_items"] == 1
        assert s["n_datasets"] == 1
        assert s["n_patients"] == 1
        assert s["n_samples"] == 1

    def test_status_active_projects(self, reg):
        reg.add_project("proj_a")
        reg.add_project("proj_b")
        s = reg.status()
        assert "proj_a" in s["active_projects"]
        assert "proj_b" in s["active_projects"]

    def test_status_phase_summary(self, reg):
        proj, ds = _setup_project_with_dataset(reg)
        reg.upsert_phase(proj, ds, "data_processing", "qc_filter", status="ready")
        reg.upsert_phase(proj, ds, "data_processing", "preprocess", status="in_progress")

        s = reg.status()
        assert proj in s["phase_summary"]
        summary = s["phase_summary"][proj]
        assert summary.get("ready", 0) == 1
        assert summary.get("in_progress", 0) == 1

    def test_status_empty_registry(self, reg):
        s = reg.status()
        assert s["n_projects"] == 0
        assert s["n_data_sources"] == 0
        assert s["n_inventory_items"] == 0
        assert s["n_datasets"] == 0
        assert s["n_patients"] == 0
        assert s["n_samples"] == 0
        assert s["active_projects"] == []
        assert s["phase_summary"] == {}


# ---------------------------------------------------------------------------
# TestEdgeCases
# ---------------------------------------------------------------------------


class TestEdgeCases:
    def test_same_inventory_item_in_multiple_datasets(self, reg):
        _make_inventory_item(reg, name="shared_item")
        reg.create_dataset("ds_a")
        reg.create_dataset("ds_b")
        reg.add_dataset_member("ds_a", "shared_item", "rna")
        reg.add_dataset_member("ds_b", "shared_item", "rna")
        members_a = reg.get_dataset_members("ds_a")
        members_b = reg.get_dataset_members("ds_b")
        assert len(members_a) == 1
        assert len(members_b) == 1
        assert members_a[0]["inventory_name"] == "shared_item"
        assert members_b[0]["inventory_name"] == "shared_item"

    def test_same_dataset_linked_to_multiple_projects(self, reg):
        reg.add_project("proj_a")
        reg.add_project("proj_b")
        reg.create_dataset("shared_ds")
        reg.link_project_dataset("proj_a", "shared_ds")
        reg.link_project_dataset("proj_b", "shared_ds")
        ds_a = reg.list_project_datasets("proj_a")
        ds_b = reg.list_project_datasets("proj_b")
        assert len(ds_a) == 1
        assert len(ds_b) == 1
        assert ds_a[0]["name"] == "shared_ds"
        assert ds_b[0]["name"] == "shared_ds"

    def test_empty_dataset_no_members(self, reg):
        reg.create_dataset("empty_ds")
        members = reg.get_dataset_members("empty_ds")
        assert members == []

    def test_delete_project_cascades_phases(self, reg):
        proj, ds = _setup_project_with_dataset(reg)
        reg.upsert_phase(proj, ds, "data_processing", "qc_filter")
        reg.delete_project(proj)
        # Phases should be gone (CASCADE)
        assert reg.get_phase(proj, "data_processing", "qc_filter") is None

    def test_delete_project_cascades_dataset_link(self, reg):
        proj, ds = _setup_project_with_dataset(reg)
        reg.delete_project(proj)
        # Dataset still exists (project_datasets CASCADE only removes the link)
        assert reg.get_dataset(ds) is not None

    def test_dataset_versioning_preserves_old_members(self, reg):
        """Members of old version are preserved independently of new version."""
        reg.create_dataset("ds_ver")
        _make_inventory_item(reg, name="item_a")
        _make_inventory_item(reg, name="item_b", uri="/data/b.h5ad", modality="protein")
        reg.add_dataset_member("ds_ver", "item_a", "rna")
        reg.add_dataset_member("ds_ver", "item_b", "protein")

        reg.bump_dataset_version("ds_ver")

        # Verify old version members are still accessible
        v1_members = reg.get_dataset_members("ds_ver", version=1)
        assert len(v1_members) == 2

        # Add new member only to v2
        _make_inventory_item(reg, name="item_c", uri="/data/c.h5ad", modality="spatial")
        reg.add_dataset_member("ds_ver", "item_c", "spatial")

        v2_members = reg.get_dataset_members("ds_ver")
        assert len(v2_members) == 3

        # v1 still has only 2
        v1_members = reg.get_dataset_members("ds_ver", version=1)
        assert len(v1_members) == 2


# ---------------------------------------------------------------------------
# TestTransactionRollback
# ---------------------------------------------------------------------------


class TestTransactionRollback:
    def test_register_inventory_item_invalid_source_no_orphan(self, reg):
        """If data_source_name is invalid, no inventory item should remain."""
        with pytest.raises(ValueError, match="DataSource.*not found"):
            reg.register_inventory_item(
                "orphan_item",
                uri="/data/orphan.h5ad",
                modality="rna",
                data_source_name="nonexistent_source",
            )
        assert reg.get_inventory_item("orphan_item") is None

    def test_add_member_invalid_inventory_no_partial(self, reg):
        """If inventory item does not exist, no dataset_member row should remain."""
        reg.create_dataset("ds_partial")
        with pytest.raises(ValueError, match="Inventory item.*not found"):
            reg.add_dataset_member("ds_partial", "nonexistent_item", "rna")
        members = reg.get_dataset_members("ds_partial")
        assert len(members) == 0

    def test_link_project_dataset_invalid_dataset_no_partial(self, reg):
        """If dataset does not exist, no project_datasets row should remain."""
        reg.add_project("proj_partial")
        with pytest.raises(ValueError, match="Dataset.*not found"):
            reg.link_project_dataset("proj_partial", "nonexistent_ds")
        results = reg.list_project_datasets("proj_partial")
        assert len(results) == 0

    def test_upsert_phase_invalid_dataset_no_partial(self, reg):
        """If dataset does not exist, no project_phases row should remain."""
        reg.add_project("proj_partial")
        with pytest.raises(ValueError, match="Dataset.*not found"):
            reg.upsert_phase("proj_partial", "nonexistent_ds", "data_processing", "qc_filter")
        phases = reg.list_phases("proj_partial")
        assert len(phases) == 0
