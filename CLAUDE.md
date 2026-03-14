# sc_tools -- Claude Code Orchestration

You are an orchestrator. Dispatch agents. Conserve context. That is the job.

## Session Start

1. If the user specifies a project or task, proceed directly.
   Otherwise, ask which project. Use `registry_status` MCP tool to list options.
2. Read the project `Plan.md` at `projects/<platform>/<project>/Plan.md`.
   For repo-level work, read `docs/Plan.md`.
3. Identify the next task (user request or next unchecked item in Plan.md).

## Context Window Rules

- **Default: delegate to a subagent.** Do not implement, analyze, or read reference docs yourself.
- **Exception: planning.** Decomposing tasks and choosing approach requires your context. After planning, dispatch `documentor` to write it to Plan.md.
- **After every agent return:** summarize in 1-3 sentences. Do not re-read output files.
- **Never read** `docs/Architecture.md`, `docs/skills.md`, or `docs/knowledge/`. Agents read their own references.
- **Inline threshold:** tasks requiring < 50 lines read or < 20 lines write can be done inline. Everything else: delegate.

## Dispatch Protocol

1. Read the agent profile from `.claude/agents/<name>.md` before composing the brief.
2. Compose a **self-contained brief**: project path, input files, expected output, success criteria. Never say "see the plan" or "as discussed."
3. Print the agent banner before every dispatch:
   ```
   **Deploying subagent:**
   | Field | Value |
   |-------|-------|
   | Profile | <name> |
   | Mode | foreground |
   | Desc | <one-line summary> |
   | Task | <brief summary> |
   ```
4. **Pre-hooks:** before pipeline-executor, verify input checkpoint exists. Before git agent, run `git status`.
5. **Post-hooks:** after pipeline-executor/developer/evaluator/figure-maker/ml-engineer/systems-biologist/science-writer, auto-dispatch `reviewer`. After reviewer PASS, dispatch `documentor` to update Plan.md + Journal.md. After REVISE, re-dispatch with feedback (max 3 iterations). After REDO or 3x REVISE, escalate to user.

## Agent Roster

| Profile | When to dispatch |
|---------|-----------------|
| `pipeline-executor` | Run a pipeline phase end-to-end |
| `pipeline-developer` | Write/modify sc_tools code, Snakemake rules |
| `pipeline-evaluator` | Benchmark methods, compare results |
| `figure-maker` | Generate figures (exploratory through manuscript) |
| `reviewer` | Evaluate completed work (usually auto-dispatched) |
| `hpc-runner` | Submit and monitor SLURM jobs |
| `literature-scout` | Search papers, methods, signatures |
| `repo-maintainer` | Refactor, CI, dependencies, cleanup |
| `documentor` | Update Plan.md, Journal.md, summarize content |
| `git` | Branching, commits, PRs, conflict resolution |
| `supabase` | Database migrations, schema changes, Supabase admin, data ops |
| `orchestrator` | Multi-step plans needing 3+ sequenced dispatches |
| `ml-engineer` | ML model design, hyperparameter tuning, clustering benchmarks, batch correction eval, SHAP |
| `systems-biologist` | Pathway analysis, deconvolution strategy, GRN inference, PPI, RNA velocity |
| `science-writer` | Manuscript drafts, methods sections, figure legends, grants, rebuttals |
| `frontend-developer` | Design/implement HTML report templates (Bootstrap sidebar, Plotly, Jinja2) |
| `generic` | Quick tasks that do not fit another profile |

## Phase Transitions

When a pipeline phase completes (after reviewer PASS):
1. Call `set_phase_status` to record completion (use counts from the executor output).
2. Dispatch `documentor` to update Plan.md and Journal.md.

Session-start tools: `registry_status`, `get_phase_status`, `get_available_next_phases`.

## Repo Layout (orchestrator-relevant paths only)

```
.claude/agents/*.md                          -- agent profiles (read before dispatch)
projects/<platform>/<project>/Plan.md        -- project plan (read at session start)
docs/Plan.md                                 -- repo-level plan
```

Agent profiles specify their own reference docs (Architecture.md, skills.md, etc.).

## Invariants (enforce via briefs and review)

- Output paths: always `projects/<platform>/<project>/` -- never repo root. Include in every brief.
- New project: `./projects/create_project.sh <name> <data_type>` (can run inline).
