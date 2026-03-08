---
name: journal-and-mission
tier: 1
description: "[Tier 1 — Core] Sync project state, record decisions, and update Mission/Journal. Use when starting a session, finishing significant work, or tracking progress in sc_tools and subprojects."
allowed-tools: Read Write Edit Glob
---

# Journal and Mission Workflow

## Files

Project documentation lives in the actual project directories. The wiki (`docs/wiki/`) symlinks to them so everything is browseable in Obsidian.

| File | Location | Role |
|------|----------|------|
| `Mission.md` | `projects/<platform>/<project>/` (symlinked into wiki) | Todo list, phase status, roadmap. **Primary work tracker.** |
| `Journal.md` | `projects/<platform>/<project>/` (symlinked into wiki) | Chronological log: dates, actions, rationale, decisions. |
| `hypotheses.md` | `docs/wiki/projects/{name}/` (wiki-native) | Scientific questions and status. Tags: `#hypothesis/investigating`, `#hypothesis/supported`. |
| `status.gen.md` | `docs/wiki/projects/{name}/` (auto-generated) | Pipeline status from registry. **Never edit manually.** |
| Root `Mission.md` | Repo root (symlinked into `docs/wiki/sc_tools/`) | Toolkit-level todo list and pipeline roadmap. |
| Root `Journal.md` | Repo root (symlinked into `docs/wiki/sc_tools/`) | Toolkit-level decision log. |

## File ownership

| Pattern | Owner | Who writes |
|---------|-------|------------|
| `*.gen.md` | System | `make wiki-sync` only |
| `*.suggest.md` | Agent | Agent proposals for human review |
| `Mission.md`, `Journal.md` | Human | Human, or agent when asked |
| `hypotheses.md` | Human | Human, or agent when asked |

## Sync (before major work)

1. Read project **`Mission.md`** — active work items and phase status
2. Read **`docs/wiki/projects/{name}/hypotheses.md`** — big picture context
3. Read **Architecture.md** — data flow, conventions, checkpoint names

## Recording (after significant work)

1. **Mission.md:** Check off completed items `[x]`, update "In Progress", add blockers/next steps
2. **Journal.md:** Append a dated entry with action, rationale, and decisions
3. **hypotheses.md:** Update status tags as hypotheses progress
4. **Phase transitions:** Use `set_phase_status` MCP tool with `n_obs`, `n_vars`, `n_samples`, then run `make wiki-sync`

## Phase transition protocol

When a phase completes:
```
set_phase_status(project, phase, "complete", n_obs=N, n_vars=N, n_samples=N, notes="...")
```
Then: `make wiki-sync` to regenerate `status.gen.md` files with updated counts.

## Suggestion workflow

- Write `.suggest.md` files (e.g. `hypotheses.suggest.md`) with proposed content
- Human promotes to `.md` (accept) or deletes (reject)

## Work mode vs plan mode

- **Work mode** (implementing/running code): Update Mission.md after significant work
- **Plan mode** (discussing/planning): Updating optional until execution starts
