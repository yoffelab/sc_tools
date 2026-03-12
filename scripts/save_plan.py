#!/usr/bin/env python3
"""Hook script: save Claude Code plan to docs/wiki/plans/ on ExitPlanMode.

Called by ~/.claude/settings.json PostToolUse hook on ExitPlanMode.
No-op when docs/wiki/plans/ does not exist in the cwd tree.

Usage (automatic via hook):
    python3 ~/.claude/hooks/save_plan.py

Usage (manual):
    python3 scripts/save_plan.py [--plan-file path/to/plan.md]
"""

from __future__ import annotations

import re
import sys
import time
from datetime import datetime
from pathlib import Path

PLANS_SOURCE = Path.home() / ".claude" / "plans"
MAX_AGE_S = 60  # seconds — look for plans written within this window


def find_repo_root() -> Path | None:
    """Walk up from cwd looking for docs/wiki/plans/."""
    cwd = Path.cwd()
    for p in [cwd, *cwd.parents]:
        if (p / "docs" / "wiki" / "plans").is_dir():
            return p
    return None


def find_recent_plan(max_age_s: int = MAX_AGE_S) -> Path | None:
    """Return most recently modified plan file written within max_age_s seconds."""
    if not PLANS_SOURCE.exists():
        return None
    now = time.time()
    candidates = [f for f in PLANS_SOURCE.glob("*.md") if now - f.stat().st_mtime <= max_age_s]
    if not candidates:
        return None
    return max(candidates, key=lambda f: f.stat().st_mtime)


def slugify(title: str) -> str:
    """Convert plan title to a URL-safe slug (max 60 chars)."""
    s = title.lower()
    s = re.sub(r"[^a-z0-9\s-]", "", s)
    s = re.sub(r"[\s-]+", "-", s.strip())
    return s[:60]


def extract_title(content: str) -> str:
    """Return text of first # heading, or 'plan' if none found."""
    for line in content.splitlines():
        line = line.strip()
        if line.startswith("# "):
            return line[2:].strip()
    return "plan"


def extract_summary(content: str) -> str:
    """Return first sentence of paragraph following ## Context heading."""
    lines = content.splitlines()
    in_context = False
    for line in lines:
        if line.strip().startswith("## Context"):
            in_context = True
            continue
        if in_context:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                break
            # First sentence of the first non-empty, non-heading line
            sent = stripped.split(".")[0]
            return (sent[:120] + "...") if len(sent) > 120 else sent
    # Fallback: first non-empty, non-heading line
    for line in lines:
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            sent = stripped.split(".")[0]
            return (sent[:120] + "...") if len(sent) > 120 else sent
    return ""


def build_output(content: str, title: str, summary: str, date: str) -> str:
    """Prepend YAML frontmatter to plan content."""
    # Escape double quotes in summary
    safe_summary = summary.replace('"', "'")
    fm = f'---\nstatus: active\ncreated: {date}\nsummary: "{safe_summary}"\n---\n'
    return fm + content


def update_mission_active_plans(
    mission_path: Path,
    plan_title: str,
    plan_slug: str,
    date: str,
) -> bool:
    """Append plan wikilink to ## Active Plans sentinel block in Mission.md.

    Returns True if the file was modified.
    """
    if not mission_path.exists():
        return False
    content = mission_path.read_text()
    if "<!-- ACTIVE_PLANS:START -->" not in content:
        return False
    filename = f"{date}-{plan_slug}"
    link = f"[[plans/untiered/{filename}|{plan_title}]]"
    entry = f"- {link} — active ({date})"
    # Skip if already referenced (idempotent)
    if filename in content:
        return False
    # Insert entry before the END sentinel
    new_content = content.replace(
        "<!-- ACTIVE_PLANS:END -->",
        f"{entry}\n<!-- ACTIVE_PLANS:END -->",
    )
    mission_path.write_text(new_content)
    return True


def main(plan_file: Path | None = None) -> int:
    """Save plan file to wiki and update Mission.md.

    Returns 0 on success or non-fatal skip; 1 on error.
    """
    if plan_file is None:
        plan_file = find_recent_plan()
    if plan_file is None:
        print("save_plan: no recent plan file found (within 60s) — skipping", file=sys.stderr)
        return 0

    repo_root = find_repo_root()
    if repo_root is None:
        print("save_plan: docs/wiki/plans/ not found in cwd tree — skipping", file=sys.stderr)
        return 0

    plans_dir = repo_root / "docs" / "wiki" / "plans" / "untiered"

    content = plan_file.read_text()
    title = extract_title(content)
    summary = extract_summary(content)
    slug = slugify(title)
    date = datetime.now().strftime("%Y-%m-%d")
    dest_name = f"{date}-{slug}.md"
    dest = plans_dir / dest_name

    output = build_output(content, title, summary, date)
    dest.write_text(output)
    print(f"save_plan: wrote {dest.relative_to(repo_root)}")

    mission = repo_root / "docs" / "Mission.md"
    changed = update_mission_active_plans(mission, title, slug, date)
    if changed:
        print("save_plan: updated Mission.md ## Active Plans")

    return 0


if __name__ == "__main__":
    # Support --plan-file override for manual/test use
    plan_file: Path | None = None
    args = sys.argv[1:]
    if "--plan-file" in args:
        idx = args.index("--plan-file")
        if idx + 1 < len(args):
            plan_file = Path(args[idx + 1])
    sys.exit(main(plan_file=plan_file))
