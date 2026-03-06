---
name: writing-plans
tier: 3
description: "[Tier 3 — Specialized] Write structured implementation plans before starting complex tasks. Use when a task spans multiple files, has architectural decisions, or needs user approval before execution. Produces a clear spec the user can approve."
allowed-tools: Read Glob Grep
---

# Writing Plans

## When to Write a Plan

Write a plan before starting when ANY of these are true:
- The task touches 3+ files
- There are multiple valid approaches (architectural decision)
- The task could be irreversible or disruptive (delete, rename, restructure)
- The user has not specified the approach and it matters
- You are not confident you fully understand the scope

Skip the plan for: single-file edits, obvious bug fixes, typo corrections.

## Plan Structure

```markdown
# Plan: <short descriptive title>

## Context
One paragraph: what problem this solves and why it is being done now.

## Approach
Why this approach was chosen over alternatives (2-3 sentences max).

## Files to Change
| Action | File | Why |
|--------|------|-----|
| Create | sc_tools/tl/new_module.py | New functionality |
| Edit | sc_tools/tl/__init__.py | Export new function |
| Edit | sc_tools/tests/test_new_module.py | Add unit tests |

## Step-by-Step
1. Write failing test for `new_function()` in `tests/test_new_module.py`
2. Implement `new_function()` in `sc_tools/tl/new_module.py`
3. Export from `sc_tools/tl/__init__.py`
4. Run `make lint` and `pytest sc_tools/tests/ -q`

## Verification
- [ ] `make lint` passes
- [ ] `pytest sc_tools/tests/ -q` — all tests pass
- [ ] New function accessible via `import sc_tools.tl`
- [ ] CI passes after push

## Out of Scope
List things that are explicitly NOT part of this plan to prevent scope creep.
```

## sc_tools Plan Checklist

Before presenting a plan:
- [ ] Read `Architecture.md` to confirm no checkpoint naming violations
- [ ] Read `Mission.md` to confirm this aligns with current priorities
- [ ] Verify the plan follows sc_tools conventions (obsm keys, checkpoint names, no project-specific code in sc_tools)
- [ ] Include a test step before the implementation step (TDD order)
- [ ] Include `make lint` in verification
- [ ] Include CI check in verification

## Exploration Before Planning

Use Read/Glob/Grep (never Edit/Write) during plan mode:
```bash
# Understand scope
Glob("sc_tools/tl/*.py")
Grep("def function_name", "sc_tools/")
Read("sc_tools/tl/__init__.py")
```

## Presenting the Plan

- State explicitly: "I have not written any code yet"
- List open questions the user must answer before execution starts
- If approach choice matters (e.g. scVI vs Harmony), present options with trade-offs
- Keep the plan in the plan file so the user can review it via `ExitPlanMode`

## Anti-Patterns

- **Planning while writing code** — finish the plan before touching files
- **Over-specifying** — the plan covers WHAT and WHY; implementation details go in code
- **Under-specifying** — steps must be concrete enough for a different agent to execute them
- **Scope creep** — if you notice a related issue, put it in "Out of Scope" and create a separate issue
