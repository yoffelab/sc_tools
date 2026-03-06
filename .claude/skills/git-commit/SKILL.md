---
name: git-commit
tier: 2
description: "[Tier 2 — Workflow] Write high-quality conventional commits: atomic scope, imperative subject, body explaining why not what, Co-Authored-By trailer. Use before every commit."
allowed-tools: Bash
---

# Git Commit

## Conventional Commit Format

```
<type>(<scope>): <subject>

<body — optional>

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
```

**Types:** `feat`, `fix`, `docs`, `style`, `refactor`, `test`, `chore`, `ci`

**Subject rules:**
- Imperative mood: "add", "fix", "remove" — not "added", "fixing", "removed"
- 72 characters max
- No period at end
- No vague words: "update", "change", "misc", "various"

**Body rules (when needed):**
- Explain WHY, not WHAT (the diff shows what)
- Wrap at 72 characters
- Separate from subject with blank line

## sc_tools Examples

```
feat(qc): add filter_spots() with modality-aware thresholds

Visium uses count-based thresholds; Xenium uses area/transcript
thresholds. Both fall back to MAD outlier detection when sample
size is small.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
```

```
fix(ingest): correct load_imc_sample() to accept processed_dir

Previously required h5ad_path which did not exist at Phase 0b time.
Now auto-discovers cells.h5ad inside processed_dir.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
```

```
ci: add GitHub Actions workflow for lint, test, and docs

Runs on push/PR to main. Jobs: lint (ruff) → test matrix
(py3.10/3.11 x ubuntu/macos) → docs (sphinx) → snakemake dry-run.
Docker build on main push only.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
```

## Pre-Commit Checklist

Before committing:
1. `make lint` passes (ruff check + format)
2. `pytest sc_tools/tests/ -q` passes (or relevant subset)
3. No files >1MB staged: `git diff --cached --stat | grep -E 'MB'`
4. No secrets or credentials in staged files
5. No `*.h5ad`, `*.ai`, `*.tiff` staged
6. Staged only specific files (not `git add .`)

## HEREDOC Template

Always use HEREDOC for multi-line messages to preserve formatting:

```bash
git commit -m "$(cat <<'EOF'
feat(scope): subject line here

Body explaining why this change was made and any important
context the diff does not make obvious.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

## Atomic Commits

One logical change per commit. If you find yourself writing "and" in the subject, split it.

**Bad:** `feat(qc): add filter_spots and fix load_imc and update tests`
**Good:** Three separate commits.

## After Committing

```bash
gh run list -L 1    # confirm CI started (sc_tools)
```
