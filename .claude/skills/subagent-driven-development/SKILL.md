---
name: subagent-driven-development
tier: 1
description: "[Tier 1 — Core] Execute an implementation plan by dispatching a fresh subagent per task with two-stage review (spec compliance, then code quality). Use for multi-task implementation plans in the current session."
---

# Subagent-Driven Development

**Core principle:** Fresh subagent per task + two-stage review (spec then quality) = high quality, fast iteration.

## When to Use

When you have an implementation plan with mostly independent tasks to execute in the current session. Contrast with dispatching-parallel-agents (for independent problems) — this is for sequential plan execution with quality gates.

## The Process (per task)

1. Dispatch **implementer subagent** with full task text and all required context
2. Answer clarifying questions **before** implementation begins
3. Implementer completes work, runs tests, self-reviews
4. Dispatch **spec compliance reviewer** subagent
   - Does the implementation match the spec? No over-building or under-building?
   - If issues found → implementer fixes → reviewer re-evaluates
5. Dispatch **code quality reviewer** subagent
   - Is the implementation well-built? Clean, tested, lint-passing?
   - If issues found → implementer addresses → reviewer re-evaluates
6. Mark task complete, move to next task
7. After all tasks: dispatch **final code reviewer** for entire implementation

## Critical Rules

- **Never skip reviews** — both gates must pass
- **Spec compliance review before code quality** — wrong order invalidates the process
- **Never proceed with unfixed issues** — the reviewer must re-evaluate after fixes
- **Same implementer subagent fixes issues** — not a new agent
- Requires git worktrees (isolated workspace); incorporates TDD practices

## Key Advantages

- No file reading overhead: controller provides full context to each subagent
- Questions surfaced before work begins, not after
- Two quality gates: spec compliance prevents over/under-building; code quality ensures clean implementation
- Context stays fresh: each subagent starts clean

## sc_tools-specific usage

Typical plan structure for a new sc_tools module:
```
Task 1: Write failing tests for new function
Task 2: Implement function to pass tests
Task 3: Add to __init__.py exports and update docs
Task 4: Update Snakemake rules if pipeline-facing
```

Each task dispatches a fresh implementer + 2 reviewers before proceeding.
