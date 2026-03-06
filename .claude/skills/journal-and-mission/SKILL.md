---
name: journal-and-mission
tier: 1
description: "[Tier 1 — Core] Sync project state, record decisions, and update Mission.md as todo list. Use when starting a session, finishing significant work, or tracking progress in sc_tools and subprojects."
allowed-tools: Read Write Edit Glob
---

# Journal and Mission Workflow

## Files

| File | Role |
|------|------|
| `Journal.md` | Full decision log: dates, actions, rationale. Can grow long. |
| `journal_summary.md` | Condensed summary of Journal.md. Kept short for context window. |
| `Mission.md` | Todo list and roadmap: checkboxes, active tasks, blockers. |

Scope: repo root + each project under `projects/<platform>/<project>/`.

## Sync (before major work)

1. Read **Mission.md** — current task and status
2. Read **journal_summary.md** (and Journal.md if more detail needed)
3. Read **Architecture.md** — data flow, conventions, checkpoint names

## Recording (after significant work)

1. **Journal.md:** Append a dated entry (action, rationale, key decisions)
2. **journal_summary.md:** Update the short summary to reflect the latest phase
3. **Mission.md:** Mark completed steps `[x]`, update active tasks, add blockers

## Work mode vs plan mode

- **Work mode** (implementing/running code): Update Mission.md after EACH prompt — no exceptions
- **Plan mode** (discussing/planning): Updating Mission optional until execution starts

## journal_summary.md format

Keep short (a few bullets or one paragraph per section):
- **Scope:** Current phase and overall goal
- **Recent:** Last major decisions and outcomes
- **Conventions:** Any key naming/path conventions relevant now
