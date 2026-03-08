"""Tests for scripts/sync_wiki.py generation functions."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Add repo root to path so we can import sync_wiki
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))


@pytest.fixture()
def phases():
    from sc_tools.pipeline import STANDARD_PHASES

    return STANDARD_PHASES


class TestGeneratePhasesMd:
    def test_generates_markdown_table(self, phases):
        from sync_wiki import generate_phases_md

        md = generate_phases_md(phases)
        assert "# Pipeline Phases" in md
        assert "| Slug |" in md
        assert "`ingest_raw`" in md
        assert "`celltype_manual`" in md

    def test_all_phases_present(self, phases):
        from sync_wiki import generate_phases_md

        md = generate_phases_md(phases)
        for slug in phases:
            assert f"`{slug}`" in md

    def test_old_codes_present(self, phases):
        from sync_wiki import generate_phases_md

        md = generate_phases_md(phases)
        assert "p0a" in md
        assert "p3.5b" in md
        assert "p6/p7" in md


class TestGenerateCheckpointsMd:
    def test_generates_sections(self, phases):
        from sync_wiki import generate_checkpoints_md

        md = generate_checkpoints_md(phases)
        assert "# Checkpoint Contracts" in md
        assert "## qc_filter" in md
        assert "**Required obs:**" in md
        assert "`sample`" in md

    def test_no_checkpoint_shown(self, phases):
        from sync_wiki import generate_checkpoints_md

        md = generate_checkpoints_md(phases)
        assert "*(no checkpoint)*" in md  # ingest_raw has no checkpoint


class TestBuildPhaseTable:
    def test_sentinel_table_format(self, phases):
        from sync_wiki import _build_phase_table

        table = _build_phase_table(phases)
        assert "| Slug | Old code |" in table
        assert "| `ingest_raw` |" in table
        lines = table.strip().split("\n")
        # Header + separator + 10 phases
        assert len(lines) == 12


class TestUpdateSentinelFile:
    def test_replaces_sentinel_content(self, tmp_path):
        from sync_wiki import update_sentinel_file

        f = tmp_path / "test.md"
        f.write_text(
            "before\n<!-- PHASE_TABLE:START -->\nold content\n<!-- PHASE_TABLE:END -->\nafter"
        )
        changed = update_sentinel_file(f, "new table")
        assert changed is True
        content = f.read_text()
        assert "new table" in content
        assert "old content" not in content
        assert "before" in content
        assert "after" in content

    def test_no_sentinel_returns_false(self, tmp_path):
        from sync_wiki import update_sentinel_file

        f = tmp_path / "test.md"
        f.write_text("no sentinels here")
        changed = update_sentinel_file(f, "new table")
        assert changed is False

    def test_missing_file_returns_false(self, tmp_path):
        from sync_wiki import update_sentinel_file

        f = tmp_path / "nonexistent.md"
        changed = update_sentinel_file(f, "new table")
        assert changed is False

    def test_idempotent(self, tmp_path):
        from sync_wiki import update_sentinel_file

        f = tmp_path / "test.md"
        f.write_text("before\n<!-- PHASE_TABLE:START -->\nstable\n<!-- PHASE_TABLE:END -->\nafter")
        changed = update_sentinel_file(f, "stable")
        assert changed is False


class TestChangelogMd:
    def test_generates_changelog(self):
        from sync_wiki import generate_changelog_md

        md = generate_changelog_md()
        assert "# Phase Changelog" in md
        # Should produce a table, a fallback message, or a registry-unavailable note
        assert "| Date |" in md or "No phase transitions" in md or "Registry not available" in md


class TestProjectStatusMd:
    def test_generates_without_registry(self):
        from sync_wiki import generate_project_status_md

        md = generate_project_status_md("test_project")
        assert "test_project" in md
        assert "Registry not available" in md or "No phase data" in md
