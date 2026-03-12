"""Unit tests for BioData-related MCP tools in sc_tools.mcp.registry_server."""

from __future__ import annotations

import pytest

mcp = pytest.importorskip("mcp", reason="mcp not installed")
sqlalchemy = pytest.importorskip("sqlalchemy", reason="sqlalchemy not installed")


@pytest.fixture()
def _patch_registry(tmp_path, monkeypatch):
    """Point the MCP server at a temp DB."""
    db_url = f"sqlite:///{tmp_path / 'test_mcp.db'}"
    monkeypatch.setenv("SC_TOOLS_REGISTRY_URL", db_url)
    # Seed a project
    from sc_tools.registry import Registry

    reg = Registry(db_url=db_url)
    reg.add_project("test_proj", platform="visium")
    return reg


@pytest.fixture()
def server(_patch_registry):
    """Import the MCP module (which creates tools at import time)."""
    from sc_tools.mcp import registry_server

    return registry_server


class TestAddSubject:
    def test_add_subject(self, server):
        result = server.add_subject("PT001", organism="human", sex="M", diagnosis="DLBCL")
        assert "PT001" in result
        assert "registered" in result.lower() or "id=" in result

    def test_add_subject_idempotent(self, server):
        server.add_subject("PT002", organism="human")
        result = server.add_subject("PT002", organism="human")
        assert "PT002" in result


class TestListSubjects:
    def test_list_empty(self, server):
        result = server.list_subjects()
        assert "No subjects" in result or "0" in result

    def test_list_after_add(self, server):
        server.add_subject("PT010", diagnosis="UC")
        result = server.list_subjects()
        assert "PT010" in result


class TestAddSample:
    def test_add_sample(self, server):
        server.add_subject("PT001")
        result = server.add_sample(
            sample_id="S001",
            subject_id="PT001",
            project_name="test_proj",
            tissue="colon",
        )
        assert "S001" in result

    def test_add_sample_no_subject(self, server):
        result = server.add_sample(
            sample_id="S002",
            project_name="test_proj",
        )
        assert "S002" in result


class TestListSamples:
    def test_list_samples(self, server):
        server.add_subject("PT001")
        server.add_sample(sample_id="S001", subject_id="PT001", project_name="test_proj")
        result = server.list_samples(project_name="test_proj")
        assert "S001" in result


class TestRegisterBiodata:
    def test_register_spatial_seq(self, server):
        result = server.register_biodata(
            project_name="test_proj",
            category="spatial_seq",
            platform="visium",
            uri="/results/adata.filtered.h5ad",
            fmt="h5ad",
        )
        assert "id=" in result

    def test_register_image(self, server):
        result = server.register_biodata(
            project_name="test_proj",
            category="image",
            platform="imc",
            uri="/data/sample1.tiff",
            fmt="tiff",
        )
        assert "id=" in result


class TestListBiodata:
    def test_list_empty(self, server):
        result = server.list_biodata(project_name="test_proj")
        assert "No BioData" in result or "0" in result

    def test_list_after_register(self, server):
        server.register_biodata(
            project_name="test_proj",
            category="spatial_seq",
            platform="xenium",
            uri="/data/xenium.h5ad",
        )
        result = server.list_biodata(project_name="test_proj")
        assert "xenium" in result


class TestProjectDataSummary:
    def test_empty_project(self, server):
        result = server.project_data_summary("test_proj")
        assert "test_proj" in result
        assert "0" in result

    def test_with_data(self, server):
        server.register_biodata(
            project_name="test_proj",
            category="spatial_seq",
            platform="visium",
            uri="/a.h5ad",
        )
        server.register_biodata(
            project_name="test_proj",
            category="image",
            platform="imc",
            uri="/b.tiff",
        )
        result = server.project_data_summary("test_proj")
        assert "2" in result  # total


class TestRegistryStatusIncludesBiodata:
    def test_status_has_new_counts(self, server):
        result = server.registry_status()
        assert "Subjects" in result
        assert "Samples" in result
        assert "BioData" in result
