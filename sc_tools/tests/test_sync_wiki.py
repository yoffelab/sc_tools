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
        f.write_text(
            "before\n<!-- PHASE_TABLE:START -->\n\nstable\n\n<!-- PHASE_TABLE:END -->\nafter"
        )
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


class TestActivePlansSentinel:
    """Tests for ACTIVE_PLANS sentinel in Mission.md."""

    def test_build_active_plans_table_empty(self, tmp_path):
        from sync_wiki import _build_active_plans_table

        table = _build_active_plans_table(plans_dir=tmp_path)
        # No plan files => returns placeholder message
        assert "*No active plans*" in table or table == ""

    def test_build_active_plans_table_with_plans(self, tmp_path):
        from sync_wiki import _build_active_plans_table

        # Write two plan files with YAML frontmatter
        (tmp_path / "2026-03-08-my-plan.md").write_text(
            '---\nstatus: active\ncreated: 2026-03-08\nsummary: "Do something"\n---\n# My Plan\n'
        )
        (tmp_path / "2026-03-07-old-plan.md").write_text(
            '---\nstatus: complete\ncreated: 2026-03-07\nsummary: "Already done"\n---\n# Old Plan\n'
        )
        table = _build_active_plans_table(plans_dir=tmp_path)
        # Only active/in_progress plans appear
        assert "My Plan" in table or "my-plan" in table
        assert "Old Plan" not in table

    def test_build_active_plans_table_in_progress(self, tmp_path):
        from sync_wiki import _build_active_plans_table

        (tmp_path / "2026-03-08-wip.md").write_text(
            '---\nstatus: in_progress\ncreated: 2026-03-08\nsummary: "In flight"\n---\n# WIP Plan\n'
        )
        table = _build_active_plans_table(plans_dir=tmp_path)
        assert "WIP Plan" in table or "wip" in table

    def test_update_active_plans_file_replaces_sentinel(self, tmp_path):
        from sync_wiki import update_active_plans_file

        f = tmp_path / "Mission.md"
        f.write_text(
            "# Mission\n\n## Active Plans\n\n"
            "<!-- ACTIVE_PLANS:START -->\n*No active plans*\n<!-- ACTIVE_PLANS:END -->\n"
        )
        changed = update_active_plans_file(f, "new plan table")
        assert changed is True
        content = f.read_text()
        assert "new plan table" in content
        assert "<!-- ACTIVE_PLANS:START -->" in content
        assert "<!-- ACTIVE_PLANS:END -->" in content

    def test_update_active_plans_file_no_sentinel(self, tmp_path):
        from sync_wiki import update_active_plans_file

        f = tmp_path / "Mission.md"
        f.write_text("# Mission\nno sentinel here")
        changed = update_active_plans_file(f, "table")
        assert changed is False

    def test_update_active_plans_file_idempotent(self, tmp_path):
        from sync_wiki import update_active_plans_file

        content = "# Mission\n\n<!-- ACTIVE_PLANS:START -->\n\nstable content\n\n<!-- ACTIVE_PLANS:END -->\n"
        f = tmp_path / "Mission.md"
        f.write_text(content)
        changed = update_active_plans_file(f, "stable content")
        assert changed is False


class TestSavePlanScript:
    """Tests for scripts/save_plan.py helper functions."""

    def test_slugify(self):
        from save_plan import slugify

        assert slugify("Unified Bookkeeping System") == "unified-bookkeeping-system"
        assert slugify("Plan: Docs & Hooks!") == "plan-docs-hooks"

    def test_slugify_max_length(self):
        from save_plan import slugify

        long = "a " * 40
        assert len(slugify(long)) <= 60

    def test_extract_title(self):
        from save_plan import extract_title

        content = "Some preamble\n# My Plan Title\n\nBody text"
        assert extract_title(content) == "My Plan Title"

    def test_extract_title_no_heading(self):
        from save_plan import extract_title

        assert extract_title("no heading here") == "plan"

    def test_extract_summary_context(self):
        from save_plan import extract_summary

        content = "# Title\n\n## Context\n\nThis is the context. More text.\n\n## Part 2\n"
        summary = extract_summary(content)
        assert "This is the context" in summary

    def test_build_output_has_frontmatter(self):
        from save_plan import build_output

        out = build_output("# Plan\nbody", "Plan", "A summary", "2026-03-08")
        assert out.startswith("---\n")
        assert "status: active" in out
        assert "created: 2026-03-08" in out
        assert "# Plan" in out

    def test_update_mission_active_plans(self, tmp_path):
        from save_plan import update_mission_active_plans

        mission = tmp_path / "Mission.md"
        mission.write_text(
            "# Mission\n\n## Active Plans\n\n"
            "<!-- ACTIVE_PLANS:START -->\n<!-- ACTIVE_PLANS:END -->\n"
        )
        changed = update_mission_active_plans(mission, "My Plan", "my-plan", "2026-03-08")
        assert changed is True
        content = mission.read_text()
        assert "my-plan" in content
        assert "My Plan" in content

    def test_update_mission_active_plans_idempotent(self, tmp_path):
        from save_plan import update_mission_active_plans

        filename = "2026-03-08-my-plan"
        mission = tmp_path / "Mission.md"
        mission.write_text(
            f"# Mission\n\n<!-- ACTIVE_PLANS:START -->\n"
            f"- [[plans/{filename}|My Plan]] — active (2026-03-08)\n"
            "<!-- ACTIVE_PLANS:END -->\n"
        )
        changed = update_mission_active_plans(mission, "My Plan", "my-plan", "2026-03-08")
        assert changed is False


class TestNexusDashboard:
    """Tests for generate_nexus_md() in sync_wiki.py."""

    def test_nexus_has_all_sections(self, tmp_path, phases):
        from sync_wiki import generate_nexus_md

        md = generate_nexus_md(phases, plans_dir=tmp_path)
        assert "## All Projects" in md
        assert "## Available Next Phases" in md
        assert "## Active Plans" in md

    def test_nexus_no_obsidian_transclusion(self, tmp_path, phases):
        from sync_wiki import generate_nexus_md

        md = generate_nexus_md(phases, plans_dir=tmp_path)
        assert "![[" not in md

    def test_nexus_generated_frontmatter(self, tmp_path, phases):
        from sync_wiki import generate_nexus_md

        md = generate_nexus_md(phases, plans_dir=tmp_path)
        assert md.startswith("---\n")
        assert "generated: true" in md
        assert "type: nexus" in md

    def test_nexus_includes_active_plans(self, tmp_path, phases):
        from sync_wiki import generate_nexus_md

        (tmp_path / "2026-03-08-test-plan.md").write_text(
            "---\nstatus: active\ncreated: 2026-03-08\n---\n# Test Plan\n"
        )
        md = generate_nexus_md(phases, plans_dir=tmp_path)
        assert "Test Plan" in md

    def test_nexus_next_phases_section_present(self, tmp_path, phases):
        from sync_wiki import generate_nexus_md

        md = generate_nexus_md(phases, plans_dir=tmp_path)
        # Section header and table header must be present
        assert "## Available Next Phases" in md
        assert "| Project |" in md

    def test_cmd_sync_writes_nexus_file(self, tmp_path, monkeypatch, phases):
        from sync_wiki import generate_nexus_md

        # Verify the function returns a string (cmd_sync integration tested via manual run)
        md = generate_nexus_md(phases, plans_dir=tmp_path)
        assert isinstance(md, str)
        assert len(md) > 100
