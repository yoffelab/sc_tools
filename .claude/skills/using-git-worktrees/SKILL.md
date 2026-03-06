---
name: using-git-worktrees
description: "[Tier 2 — Workflow] Create an isolated git worktree before starting feature work or executing an implementation plan. Use when working on multiple branches simultaneously or to protect the main working tree."
tier: 2
---

# Using Git Worktrees

**Core principle:** Systematic directory selection + safety verification = reliable isolation.

**Announce at start:** "I'm using the using-git-worktrees skill to set up an isolated workspace."

## Directory Selection Priority

1. Check for existing `.worktrees/` or `worktrees/` — use whichever exists (`.worktrees/` wins if both)
2. Check `CLAUDE.md` for a worktree directory preference
3. If neither: ask the user

## Safety Verification (MANDATORY for project-local)

Before creating a project-local worktree, verify it is gitignored:
```bash
git check-ignore -q .worktrees 2>/dev/null
```
If NOT ignored: add to `.gitignore`, commit, then proceed.

## Creation Steps

```bash
# 1. Create worktree with new branch
git worktree add .worktrees/feature-name -b feature/feature-name
cd .worktrees/feature-name

# 2. Run project setup
pip install -e ".[dev]"   # for sc_tools

# 3. Verify clean baseline
pytest sc_tools/tests/ -v --tb=short -q
```

Report failures before proceeding — never start with a broken baseline.

## Quick Reference

| Situation | Action |
|-----------|--------|
| `.worktrees/` exists | Use it (verify ignored) |
| `worktrees/` exists | Use it (verify ignored) |
| Neither exists | Check CLAUDE.md → ask user |
| Not ignored | Add to .gitignore + commit first |
| Baseline tests fail | Report failures + ask whether to proceed |

## sc_tools Notes

- `.worktrees/` is already in `.gitignore` — safe to use directly
- Run `make lint` after setup to confirm ruff is configured correctly
- For HPC work: create worktree locally, rsync to scratch before submitting SLURM jobs

## Red Flags

- Creating worktree without verifying gitignore (project-local)
- Skipping baseline test verification
- Assuming directory location without checking
