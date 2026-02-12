---
name: journal-and-mission-workflow
description: Maintain journal_summary.md as condensed long-term memory and Mission.md as the todo list; in work mode, update Mission.md after each prompt. Use when syncing project state, recording decisions, or tracking progress in sc_tools and subprojects.
---

# Journal Summary and Mission-as-Todo Workflow

## Purpose

- **Journal.md** and **journal_summary.md** together form long-term memory of what was done and why. The summary shortens context and speeds up sync.
- **Mission.md** is the single source of truth for progress: a living todo list (steps, checkboxes, current priority).
- In **work mode** (executing tasks, not just planning), the agent updates Mission.md after each prompt so progress is always current.

## Files

| File | Role |
|------|------|
| `Journal.md` | Full decision log: dates, actions, rationale. Can grow long. |
| `journal_summary.md` | Condensed summary of Journal.md (phases, key decisions, outcomes). Kept short for context. |
| `Mission.md` | Todo list and roadmap: completed steps, active tasks, blockers. |

**Scope:** Use at repo root and in every project that has a Journal.md (e.g. `projects/imc/lymph_dlbcl/`, `projects/visium/ggo_visium/`). Each has its own Mission.md and Journal.md; each should have a journal_summary.md.

## Sync loop (before major work)

1. Read **mission.md** to see current task and status.
2. Read **journal_summary.md** (and Journal.md if detail needed) for recent context.
3. Read **architecture.md** for data flow and conventions.

## Recording (after significant work)

1. **Journal.md:** Append a dated entry: action, rationale, key decisions.
2. **journal_summary.md:** Update the summary so it reflects the latest phase: a short paragraph or bullet list of what was decided and done (no need to copy every Journal entry).
3. **Mission.md:** Treat as the todo list. Mark completed steps (e.g. `[x]`), set "In progress" where relevant, add or adjust next steps and blockers.

## Work mode vs plan mode

- **Work mode:** You are implementing, running scripts, or making concrete changes. After each prompt in work mode, update **Mission.md** (check off completed items, set current step, note blockers). Do not skip this.
- **Plan mode:** You are only discussing, outlining, or planning. Updating Mission.md is optional until execution starts.

## journal_summary.md format

Keep it short (one to a few short sections). Prefer:

- **Repo/project summary:** One short paragraph on scope and current phase.
- **Recent phase:** Bullets or one paragraph on the last major decisions and outcomes (derived from Journal.md).
- Optionally: **Key conventions** (e.g. panel names, checkpoint names) in one line each.

Update the summary whenever Journal.md gets a new substantive entry; do not duplicate full log text.

## Mission.md as todo list

- Use checkboxes: `[ ]` for todo, `[x]` for done.
- Keep a clear "Current Status" or "Active Tasks" section at the top or in a dedicated section.
- After each prompt in work mode: re-read Mission, then update checkboxes and status lines so the file reflects actual progress.
