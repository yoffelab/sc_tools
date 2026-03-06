---
name: executing-plans
tier: 3
description: "[Tier 3 — Specialized] Execute an approved implementation plan step by step: one step at a time, verify after each step, never skip verification. Use after a plan has been approved by the user."
allowed-tools: Read Write Edit Bash Glob Grep
---

# Executing Plans

## Core Principle

Execute one step at a time. Verify before moving to the next. Never skip steps.

## Execution Loop

For each step in the approved plan:

1. **Announce** — "Executing step N: [description]"
2. **Execute** — make the change (write/edit file, run command)
3. **Verify** — run the relevant check (lint, test, dry-run)
4. **Record** — note what was done and whether it passed
5. **Stop if blocked** — do not attempt the next step if this one fails

```
Step 1: Write failing test → run pytest → confirm RED
Step 2: Implement function → run pytest → confirm GREEN
Step 3: Export from __init__.py → run import test → confirm accessible
Step 4: make lint → confirm clean
Step 5: git add specific files → git commit → gh run list -L 1
```

## When to Pause

Stop and check with the user when:
- A step produces unexpected output (error that was not in the plan)
- A file already contains code that conflicts with the plan
- You realize a step in the plan was underspecified and you must make a choice
- A test was supposed to fail but passed (may indicate the test is wrong)

Do NOT silently make decisions that deviate from the approved plan.

## Verification After Each Step

| Step type | Verification command |
|-----------|----------------------|
| Write new Python file | `python -c "import sc_tools.<module>"` |
| Add function | `pytest sc_tools/tests/test_<module>.py -q` |
| Edit __init__.py | `python -c "from sc_tools.tl import <function>"` |
| Edit Snakefile | `snakemake --snakefile <path> --dry-run all` |
| Any code change | `make lint` |
| Before commit | `pytest sc_tools/tests/ -q && make lint` |
| After push | `gh run list -L 1` |

## sc_tools Execution Order

Always follow this order when adding a new sc_tools function:

1. Write failing test in `sc_tools/tests/test_<module>.py`
2. Run test — confirm it fails (RED)
3. Implement function in `sc_tools/<module>/<file>.py`
4. Run test — confirm it passes (GREEN)
5. Export from `sc_tools/<module>/__init__.py`
6. `make lint` — confirm clean
7. Stage specific files with `git add`
8. Commit with conventional commit message
9. Push and `gh run list -L 1`

## Anti-Patterns

- **Batch execution** — running all steps without intermediate verification leads to cascading failures that are hard to untangle
- **Skipping the RED step in TDD** — if you do not confirm the test fails first, you may be testing the wrong thing
- **Continuing after a failed verification** — fix the failure before proceeding; do not skip ahead
- **Adding files outside the plan scope** — stick to the files listed in the plan; note any necessary additions for user review
- **Silently fixing "minor" issues** — if you fix something not in the plan, announce it

## Post-Execution Report

After completing all steps:

```
Execution complete.

Steps completed:
- [x] Step 1: Written test — PASSED
- [x] Step 2: Implemented function — PASSED
- [x] Step 3: Updated __init__.py — PASSED
- [x] Step 4: make lint — CLEAN
- [x] Step 5: Committed and pushed — CI started (run #XXXX)

Deviations from plan:
- None  (or: list any)

CI status: in progress — run `gh run list -L 1` to check
```
