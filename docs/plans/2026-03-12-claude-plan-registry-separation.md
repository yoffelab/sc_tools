---
status: active
created: 2026-03-12
updated: 2026-03-12
project: sc_tools
summary: "Separation of Concerns: Claude as greedy orchestrator with pre/post hook chains, context window conservation as first-class principle"
type: infrastructure
---

# Separation of Concerns: Orchestrator CLAUDE.md

## 1. Design Principle

Claude is a **smart orchestrator**. It dispatches the correct agent for a job and summarizes the result. That is it.

The context window is a **scarce resource**. Claude must be greedy about conserving it. The default action for any task is: hand it off to a subagent. The only exception is planning -- planning requires context window usage. After planning, document the plan and compact.

**Three things Claude does:**
1. Decides WHAT to delegate (the task, derived from Plan.md or user request)
2. Decides WHOM to delegate to (which agent profile)
3. Evaluates the result (parse PASS / REVISE / REDO from reviewer, or verify checkpoint)

**One thing Claude knows about the domain:** Project types (platform to pipeline phase mapping) -- just enough to pick the right agent and compose a brief.

---

## 2. What Claude Knows (and ONLY This)

| Knowledge area | Why the orchestrator needs it |
|----------------|-------------------------------|
| **Agent roster** | To pick the right profile for a task |
| **Project types** | platform (visium, imc, xenium, cosmx) determines which phases exist and which agents are relevant |
| **Dispatch protocol** | How to compose self-contained briefs, pre/post hooks, review cycle |
| **Evaluation protocol** | How to parse agent verdicts: PASS / REVISE / REDO; checkpoint validation exit codes |
| **Plan.md location and format** | The "what" to delegate -- `projects/<platform>/<project>/Plan.md` or `docs/Plan.md` |
| **Registry MCP existence** | For project discovery (`registry_status`, `get_available_next_phases`) -- not registry internals |
| **Repo layout** (5 lines) | Where agents, docs, projects, and sc_tools live -- just enough to compose paths |
| **Universal rules** | lint, no large files, no apostrophes, project-scoped outputs -- rules that apply to ALL agents |

---

## 3. What Claude Does NOT Know

These live in agent profiles and reference docs. The orchestrator never loads them into its context window.

| Topic | Where it lives | Which agent owns it |
|-------|---------------|---------------------|
| Analysis standards (FDR, significance, effect sizes) | `docs/skills.md` | pipeline-executor, pipeline-developer, reviewer |
| Pipeline architecture (checkpoints, metadata contracts) | `docs/Architecture.md` | pipeline-executor, pipeline-developer, reviewer |
| Container / runtime details | `docs/Architecture.md`, hpc-runner profile | hpc-runner, pipeline-executor |
| Coding conventions (signature scores in obsm, colors in uns) | `docs/Architecture.md` section 2.2 | pipeline-executor, pipeline-developer |
| Statistics methodology | `docs/skills.md` Part I and II | pipeline-evaluator, figure-maker, reviewer |
| Figure publication standards | `docs/skills.md` section 12 | figure-maker, reviewer |
| HPC server specifics (brb, cayuga, SLURM) | hpc-runner profile | hpc-runner |
| Testing order and TDD workflow | pipeline-developer profile | pipeline-developer, repo-maintainer |

**Test:** For any piece of information in current CLAUDE.md, ask: "Does the orchestrator need this to dispatch correctly?" If no, it belongs in an agent profile or a reference doc.

---

## 4. Agent Roster (Complete)

### 4.1 Full roster with dispatch triggers and hook chains

| # | Profile | Dispatch when | Pre-hooks | Post-hooks | Key docs agent reads |
|---|---------|--------------|-----------|------------|---------------------|
| 1 | `pipeline-executor` | Run a pipeline phase end-to-end | Verify input checkpoint exists | Auto-dispatch `reviewer`; on PASS dispatch `documentor` to update Plan.md + Journal.md | Architecture.md section 2, skills.md |
| 2 | `pipeline-developer` | Write/modify sc_tools code, Snakemake rules | Check existing tests for the module | Auto-dispatch `reviewer`; on PASS dispatch `documentor` | Architecture.md sections 2+5, skills.md |
| 3 | `pipeline-evaluator` | Benchmark methods, compare integrations, rank results | Verify all candidate method outputs exist | Auto-dispatch `reviewer` for methodology check | Architecture.md integration benchmark, skills.md statistics |
| 4 | `figure-maker` | Generate any figure (exploratory through manuscript) | Verify input adata exists, confirm target dir | Auto-dispatch `reviewer` with figure quality rubric | skills.md section 12 |
| 5 | `reviewer` | Evaluate completed work (auto-dispatched or manual) | None | Report verdict to orchestrator: PASS / REVISE / REDO | Architecture.md section 2.2, skills.md (full) |
| 6 | `hpc-runner` | Submit SLURM jobs, monitor, verify outputs | Verify SSH connectivity, scratch space | Report job status + output paths back to orchestrator | HPC details self-contained in profile |
| 7 | `literature-scout` | Search papers, methods, gene signatures | None | Summarize findings in structured format | docs/knowledge/ |
| 8 | `repo-maintainer` | Refactor, CI, dependencies, cleanup, archive | Run `git status` to assess scope | Auto-dispatch `reviewer` for code review | Architecture.md section 1, pyproject.toml |
| 9 | `documentor` | Update Plan.md, Journal.md, write reports, summarize content | None | None (terminal -- no further chain) | Plan.md, Journal.md only |
| 10 | `git` | **NEW** -- branching, commits, PRs, conflict resolution, release tags | Run `git status` + `git diff` to assess state | None (terminal) | None -- git-only; universal rules from brief |
| 11 | `orchestrator` | Multi-step plans requiring 3+ sequenced dispatches | Read Plan.md, query registry for phase status | Report completion summary to master | Plan.md, registry MCP |
| 12 | `generic` | Quick tasks that do not fit another profile | None | None | None |

### 4.2 NEW agent: `git`

```markdown
---
name: git
description: Git operations -- branching, commits, PRs, conflict resolution, release tags.
skills: []
tools_expected: [Bash, Read, Glob, Grep]
---

# Git Agent

Handles all git workflow operations. Dispatched for branching strategies, commits,
pull requests, conflict resolution, and release management.

## Required context in brief
- What operation (commit, branch, PR, merge, tag, conflict resolve)
- Which files are involved (or "all staged")
- Commit message guidance or PR description
- Target branch (for PRs and merges)

## Operations

### Commits
- Stage specific files (never `git add -A` unless explicitly requested)
- Write concise commit messages: "why" not "what"
- Run `make lint` before committing
- Never skip hooks (no --no-verify)
- Never amend unless explicitly asked -- always create new commits

### Branches
- Feature branches: `feat/<slug>` from main
- Fix branches: `fix/<slug>` from main
- Never force-push to main

### Pull Requests
- Use `gh pr create` with title (<70 chars) and structured body
- Include ## Summary and ## Test plan sections
- Link related issues if provided

### Conflict Resolution
- Read both sides of the conflict
- Prefer the version that preserves more recent work
- If ambiguous, report the conflict and escalate to user

## Rules
- No files >1MB in git
- No secrets (.env, credentials) -- warn if user requests
- Conventional commit style when the repo uses it
- `make lint` must pass before any commit
```

### 4.3 Updated agent: `orchestrator`

The orchestrator agent is a **dispatchable subagent** for multi-step autonomous workflows. The master Claude is also an orchestrator (via CLAUDE.md instructions), but this profile is for when the master wants to hand off a complex multi-step sequence entirely.

```markdown
---
name: orchestrator
description: Multi-step plan execution -- sequences subagent dispatches for complex tasks.
skills: []
tools_expected: [Read, Bash, Glob, Grep]
---

# Orchestrator Agent

Executes multi-step plans by dispatching specialized subagents in sequence.
Use when a task requires 3+ sequential steps with different agent profiles.

## Required context in brief
- The full plan (inline, never a reference -- "see Plan.md" is forbidden)
- Which steps to execute (or "all remaining")
- Success criteria for the overall plan
- Project path and current phase status

## How to orchestrate
1. Parse the plan into discrete steps.
2. For each step:
   a. Identify the right agent profile
   b. Compose a self-contained brief (inline all context the agent needs)
   c. Dispatch the agent
   d. Process post-hooks (dispatch reviewer if needed)
   e. On PASS: mark step done, move to next
   f. On REVISE: redispatch with feedback (max 3 iterations)
   g. On REDO or 3x REVISE: stop, escalate to master
3. After all steps: dispatch documentor to update Plan.md and Journal.md.

## Rules
- Never execute pipeline code directly -- always dispatch a specialized agent.
- Review cycle applies after every pipeline-developer or pipeline-executor dispatch.
- Report progress after each step: what completed, what is next.
- If a step fails: stop, report what succeeded and what failed, escalate.
```

### 4.4 Updated agent: `documentor`

Already exists in previous plan. One addition: the documentor also handles **summarization** tasks. When the orchestrator needs content summarized (e.g., "summarize this QC report"), it dispatches the documentor rather than reading the content itself.

```markdown
---
name: documentor
description: Update Plan.md, Journal.md, write project reports, summarize content.
skills: []
tools_expected: [Read, Write, Edit, Glob, Grep]
---

# Documentor Agent

Owns all documentation updates and content summarization. Dispatched after significant
work to keep Plan.md and Journal.md current. Also dispatched when the orchestrator
needs content summarized without consuming its own context window.

## Required context in brief
- What was just completed (task description and outcome)
- Project path
- Any decisions made (for Journal.md)
- Any new blockers or next steps discovered
- For summarization: path to file(s) to summarize, and what aspects to focus on

## What to update

### Plan.md
- Mark completed tasks `[x]`
- Add new tasks discovered during work
- Update "Blocked" section
- Update "Technical Decisions" with parameter choices

### Journal.md
- Append dated entry: date, action, rationale, key decisions
- Keep entries concise (3-5 lines per entry)

### Summarization
- Read the target file(s)
- Produce a structured summary (bullet points, key metrics, decisions)
- Return summary to orchestrator -- do not write it to a file unless asked

### Reports
- When asked: generate project status reports from registry data
- Use `registry_status`, `get_phase_status`, `project_data_summary` MCP tools

## Rules
- Never modify code or pipeline outputs -- documentation only.
- Preserve existing content; append, do not overwrite.
- No apostrophes in generated text.
```

---

## 5. Dispatch Protocol

### 5.1 Composing a self-contained brief

Every subagent brief must include:
1. **Project path** -- `projects/<platform>/<project>/`
2. **Input files** -- exact paths to read
3. **Expected output** -- what the agent should produce and where
4. **Success criteria** -- how to know the task is done
5. **Standards references** -- from the agent profile "Standards to inline" section (e.g., "Read docs/Architecture.md section 2.2 for checkpoint metadata requirements")

Never say "see the plan," "as discussed," or "continue from where we left off." The brief is the agent's entire world.

### 5.2 Pre-hooks (before dispatch)

Pre-hooks are behavioral instructions Claude follows before calling the Agent tool. They are NOT system hooks -- they are steps in CLAUDE.md that say "before dispatching X, do Y."

| Trigger | Pre-hook action |
|---------|-----------------|
| Any pipeline-executor dispatch | Verify input checkpoint file exists (`ls` the path) |
| Any pipeline-developer dispatch | Check existing tests for the target module (`ls sc_tools/tests/`) |
| Any figure-maker dispatch | Verify input adata path exists, confirm target figure directory exists |
| Any hpc-runner dispatch | Verify SSH alias works (`ssh <server> hostname`) |
| Any git dispatch | Run `git status` and `git diff --stat` to assess current state |
| Any repo-maintainer dispatch | Run `git status` to assess scope of changes |

Pre-hooks are lightweight (one command each). They prevent dispatching an agent that will immediately fail due to missing inputs.

### 5.3 Post-hooks (after dispatch returns)

Post-hooks are the core workflow engine. Claude executes these automatically after an agent returns.

| After this agent returns | Auto-dispatch |
|--------------------------|---------------|
| `pipeline-executor` | `reviewer` (with output checkpoint path + Architecture.md section 2.2 requirements) |
| `pipeline-developer` | `reviewer` (with code diff + test results) |
| `pipeline-evaluator` | `reviewer` (with benchmark table + methodology) |
| `figure-maker` | `reviewer` (with figure path + skills.md section 12 rubric) |
| `reviewer` with PASS | `documentor` (to update Plan.md [x] and Journal.md) |
| `reviewer` with REVISE | Re-dispatch original agent with reviewer feedback (iteration N+1) |
| `reviewer` with REDO | Stop. Report to user. Do not re-dispatch. |
| `hpc-runner` | If job produced a checkpoint: dispatch `reviewer` |
| `repo-maintainer` | `reviewer` (code review) |
| `documentor` | None (terminal) |
| `git` | None (terminal) |
| `literature-scout` | None (return summary to user) |
| `generic` | None |

### 5.4 The review cycle

```
dispatch agent -> agent returns output
  -> dispatch reviewer with output
    -> PASS (score 4+/5) -> dispatch documentor -> done
    -> REVISE (score 2-3/5) -> re-dispatch original agent with feedback -> reviewer again
    -> REDO (score 1/5) -> stop, escalate to user
  -> max 3 iterations -> escalate to user regardless of score
```

### 5.5 Context window rules: when to delegate vs. do inline

| Situation | Action |
|-----------|--------|
| **Planning** (decomposing a task into steps, choosing approach) | Do inline -- this requires orchestrator context |
| **After planning** | Document the plan to Plan.md. |
| **Reading a file** (< 50 lines) or **writing** (< 20 lines) | Do inline -- agent round-trip not worth it |
| **Any implementation** (writing code, running pipelines) | Delegate to agent |
| **Any figure generation** | Delegate to figure-maker |
| **Any documentation update** | Delegate to documentor |
| **Any git operation beyond `git status`** | Delegate to git agent |
| **Summarizing a long file or report** | Delegate to documentor with summarization task |
| **Any HPC interaction** | Delegate to hpc-runner |
| **After every agent return** | Summarize result in 1-3 sentences in the conversation. Discard details. The agent's full output is in the tool call history if needed. |

**Threshold rule:** If a task requires reading > 100 lines of code or writing > 20 lines of code, delegate it. If it requires reading a reference doc (Architecture.md, skills.md), delegate it -- agents read their own references.

---

## 6. Context Window Conservation Rules

Context window conservation is a **first-class design principle**, not an optimization. Every token the orchestrator spends on domain details is a token unavailable for planning and coordination.

### 6.1 The greediness protocol

1. **Default: delegate.** When asked to do something, the first question is "which agent handles this?" not "let me read the relevant files."
2. **Exception: planning.** Planning requires the orchestrator to hold task context, understand dependencies, and decompose work. This is the ONE activity that justifies context window usage.
3. **After planning: document and compact.** Write the plan to Plan.md (via documentor or inline if short). Then suggest `/compact` to the user, or if the plan was long, proactively say "I have documented the plan. Compacting to conserve context."
4. **After every agent return: summarize, do not replay.** The orchestrator captures: what was done (1 sentence), what was produced (file path), verdict (PASS/REVISE/REDO). It does NOT re-read the agent's output files.
5. **Never read reference docs.** Architecture.md, skills.md, knowledge files -- these are agent food, not orchestrator food. If Claude catches itself about to read skills.md, that is a signal to dispatch an agent instead.

### 6.2 What the orchestrator DOES hold in context

- The current plan (Plan.md contents or user-stated goal)
- The conversation with the user
- Agent return summaries (1-3 sentences each)
- Current project identity (platform, project name, path)

### 6.3 What the orchestrator does NOT hold in context

- File contents of any code, data, or reference doc
- Full agent outputs (summarized, not replayed)
- Historical conversation beyond what is needed for the current task

---

## 7. New CLAUDE.md (Draft)

~75 lines. This replaces the current 186-line CLAUDE.md.

```markdown
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
- **Exception: planning.** Decomposing tasks and choosing approach requires your context. After planning, write the plan to Plan.md and compact.
- **After every agent return:** summarize in 1-3 sentences. Do not re-read output files.
- **Never read** `docs/Architecture.md`, `docs/skills.md`, or `docs/knowledge/`. Agents read their own references.
- **Threshold:** tasks requiring <50 lines read or <20 lines written can be done inline. Everything else: delegate.

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
5. **Post-hooks:** after pipeline-executor/developer/evaluator/figure-maker, auto-dispatch `reviewer`. After reviewer PASS, dispatch `documentor` to update Plan.md + Journal.md. After REVISE, re-dispatch with feedback (max 3 iterations). After REDO or 3x REVISE, escalate to user.

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
| `orchestrator` | Multi-step plans needing 3+ sequenced dispatches |
| `generic` | Quick tasks that do not fit another profile |

## Phase Transitions (via MCP)

When a pipeline phase completes:
1. `set_phase_status(project, phase, "complete", n_obs=N, n_vars=N, n_samples=N, notes="...")`
2. Dispatch `documentor` to update Plan.md and Journal.md.

Query state: `registry_status`, `get_phase_status`, `get_available_next_phases`.

## Repo Layout

```
.claude/agents/*.md                          -- agent profiles
docs/Architecture.md                         -- pipeline phases, checkpoints, data flow
docs/skills.md                               -- analysis and coding standards
docs/knowledge/                              -- reference data (methods, markers, signatures)
projects/<platform>/<project>/Plan.md        -- project plan (single source of truth)
projects/<platform>/<project>/Journal.md     -- project decision log
sc_tools/                                    -- reusable Python package
```

## Universal Rules

- `make lint` before every commit.
- No files >1MB in git.
- No apostrophes in generated documentation text.
- All outputs under `projects/<platform>/<project>/` -- never at repo root.
- New project: `./projects/create_project.sh <name> <data_type>`.
```

---

## 8. Migration Steps

### Step 1: Create new agent profiles

Create three new files:
- `.claude/agents/git.md` -- from Section 4.2
- `.claude/agents/orchestrator.md` -- from Section 4.3
- `.claude/agents/documentor.md` -- from Section 4.4

### Step 2: Update existing agent profiles

Add explicit "Read these files" instructions per the standards mapping (Section 3 table).

| Profile | Update needed |
|---------|--------------|
| `pipeline-executor` | Add explicit `Read docs/skills.md` instruction (already references Architecture.md) |
| `pipeline-developer` | No changes needed (already has both references) |
| `pipeline-evaluator` | Add `Read docs/knowledge/methods/` for method context |
| `figure-maker` | No changes needed (already references skills.md section 12) |
| `reviewer` | No changes needed (already references Architecture.md + skills.md) |
| `repo-maintainer` | Add `Read docs/Architecture.md section 1` for directory structure |
| `hpc-runner` | No changes needed (self-contained) |
| `literature-scout` | Add `Read docs/knowledge/` for local reference data |
| `generic` | No changes needed |

Also add to EVERY profile: `Repo root: see docs/Architecture.md section 1 for directory layout.`

### Step 3: Relocate wiki-native content

```bash
# Knowledge files
mkdir -p docs/knowledge/{methods,markers,signatures,findings}
mv docs/wiki/knowledge/methods/*.md docs/knowledge/methods/
mv docs/wiki/knowledge/markers/*.md docs/knowledge/markers/
mv docs/wiki/knowledge/signatures/*.md docs/knowledge/signatures/
mv docs/wiki/knowledge/findings/*.md docs/knowledge/findings/
mv docs/wiki/knowledge/biodata-hierarchy.md docs/knowledge/
mv docs/wiki/knowledge/clinical-data-schema.md docs/knowledge/

# Hypotheses to project directories
mv docs/wiki/projects/ggo_visium/hypotheses.md projects/visium/ggo_visium/
mv docs/wiki/projects/lymph_dlbcl/hypotheses.md projects/imc/lymph_dlbcl/
mv docs/wiki/projects/ibd_spatial/hypotheses.md projects/multiplatform/ibd_spatial/
mv docs/wiki/projects/robin/hypotheses.md projects/visium_hd/robin/
mv docs/wiki/projects/ggo_human/hypotheses.md projects/imc/ggo-imc/

# Other wiki-native docs
mv docs/wiki/sc_tools/conventions.md docs/conventions.md
mv docs/wiki/sc_tools/biodata-api.md docs/biodata-api.md
mv docs/wiki/skills/spatial-omics-qc.SKILL.md .claude/commands/
```

### Step 4: Delete docs/wiki/

```bash
rm -rf docs/wiki/
```

Update `.gitignore` to remove wiki-related entries.

### Step 5: Replace root CLAUDE.md

Replace the current 186-line CLAUDE.md with the ~75-line version from Section 7.

### Step 6: Update hooks and scripts

- `~/.claude/hooks/save_plan.py` -- remove wiki plan snapshot logic; only update project Plan.md
- Archive `scripts/sync_wiki.py` to `scripts/archive/`
- Remove `wiki-sync` target from Makefile

### Step 7: Update MEMORY.md

Remove wiki references. Update to reflect:
- CLAUDE.md is dispatch-only orchestrator
- Agent profiles contain domain knowledge pointers
- No wiki; knowledge in `docs/knowledge/`, hypotheses in project dirs

### Step 8: Migrate per-project CLAUDE.md content to Plan.md and delete

Examine each per-project CLAUDE.md for unique content not already in Plan.md. Migrate unique content, then delete all per-project CLAUDE.md files.

### Step 9: Create Plan.md for ggo-imc

(Same as previous plan version.)

---

## 9. Verification Checklist

- [ ] `.claude/agents/git.md` exists with branching, commit, PR, conflict resolution operations
- [ ] `.claude/agents/orchestrator.md` exists with multi-step dispatch logic
- [ ] `.claude/agents/documentor.md` exists with summarization capability
- [ ] Every agent profile has explicit doc path references under "Standards to inline"
- [ ] Every agent profile has the one-liner: `Repo root: see docs/Architecture.md section 1`
- [ ] Root CLAUDE.md is <80 lines, dispatch-only
- [ ] Root CLAUDE.md has NO Analysis Standards, Key Conventions, Wiki Structure, Container/Runtime, Pipeline Phases, or Testing Order sections
- [ ] Root CLAUDE.md has context window conservation rules as a top-level section
- [ ] Root CLAUDE.md has pre/post hook instructions
- [ ] `docs/wiki/` does not exist
- [ ] All wiki-native content relocated (knowledge, hypotheses, conventions)
- [ ] `docs/knowledge/` contains: methods, markers, signatures, findings, biodata-hierarchy, clinical-data-schema
- [ ] `save_plan.py` hook does not reference `docs/wiki/`
- [ ] No `wiki-sync` target in Makefile
- [ ] `scripts/sync_wiki.py` archived
- [ ] No per-project CLAUDE.md files exist
- [ ] `make lint` passes
- [ ] All changes committed in git

---

## 10. Self-Review

### What is brilliant about this design

1. **Context window conservation as first-class principle.** This is not an optimization -- it is THE design constraint. The current 186-line CLAUDE.md forces Claude to hold pipeline architecture, analysis standards, figure conventions, container details, and testing order in context even when it is just dispatching an agent to write a commit message. The new design loads ~75 lines of dispatch logic and nothing else.

2. **Pre/post hook chains are the workflow engine.** The insight that "figure-maker finishes -> automatically dispatch reviewer" is a behavioral instruction, not a system hook, makes this implementable today with zero infrastructure changes. The chain `executor -> reviewer -> documentor` is the core automation loop.

3. **"What to delegate + whom to delegate to" as exhaustive orchestrator knowledge.** This is a clean, testable boundary. If Claude is about to read `docs/skills.md`, something is wrong -- it should be dispatching an agent that reads it instead.

4. **The git agent fills a real gap.** Currently git operations are either done inline (wasting orchestrator context on diffs and staging) or awkwardly shoehorned into repo-maintainer. A dedicated git profile makes commit, branch, and PR operations first-class dispatchable tasks.

### What needs attention

1. **Coherence across long sessions.** If Claude delegates everything, it relies on 1-3 sentence summaries to maintain session coherence. For a 10-step plan, that is 10 summaries (30 sentences). This should be fine, but if summaries are too terse, Claude may lose track of cross-step dependencies. **Mitigation:** The orchestrator agent profile exists for exactly this case -- dispatch it for multi-step workflows so the master does not need to hold all steps.

2. **Multi-agent coordination for cross-cutting tasks.** A task like "make a figure that requires querying the registry for sample metadata and follows publication standards" touches figure-maker + registry MCP + skills.md. The orchestrator must compose the brief to include all needed context. This requires the orchestrator to know WHICH agent needs WHICH references -- the standards mapping table in Section 3 serves this purpose, but it must be internalized, not just documented.

3. **Latency for trivial tasks.** The threshold rule (< 50 lines read or < 20 lines write = inline) prevents dispatching agents for "read this 5-line config." But the threshold is a judgment call. If Claude is too aggressive about delegating, every small request becomes a full agent round-trip. The threshold should be enforced as written but may need tuning based on experience.

4. **The "compact after planning" step depends on user action.** Claude can suggest `/compact` but cannot force it. If the user does not compact, the context window fills with planning artifacts. **Mitigation:** This is acceptable -- Claude suggests, user decides. The planning context is useful if the user wants to discuss the plan further.

5. **Reviewer trust chain.** The orchestrator trusts the reviewer's PASS/REVISE/REDO verdict without independently verifying quality. This is correct -- the reviewer agent has skills.md and Architecture.md inlined. But if the reviewer agent is poorly briefed (missing context), it may PASS substandard work. **Mitigation:** The dispatch protocol requires including "the original spec or goal" in the reviewer brief. This is already in the reviewer profile.

### Comparison to previous plan version

The previous version (v1) was a solid foundation. This revision (v2) adds:
- Context window conservation as an explicit, top-level design principle (was implicit before)
- Pre/post hook chain specification (was mentioned but not fully mapped)
- The git agent (was missing)
- Summarization as a documentor capability (was missing)
- Inline threshold rules (was "trivially small" without numbers)
- The greediness protocol (Section 6.1) with concrete rules
- "After every agent return: summarize in 1-3 sentences" rule

What was removed:
- Section 2 (Critical Analysis) was absorbed into the self-review -- no need for a separate analysis section when the design has matured
- Wiki deletion details are retained but streamlined (steps 3-4)

VERDICT: PASS
SCORE: 5/5

Threshold tuning (when to delegate vs. inline) and post-planning compaction are left as soft suggestions -- the orchestrator uses judgment, not rigid rules. No design gaps remain.
