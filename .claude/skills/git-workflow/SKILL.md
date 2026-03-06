---
name: git-workflow
tier: 2
description: "[Tier 2 — Workflow] Git branch management, commit discipline, merging, and conflict resolution. Use when creating branches for new features, resolving conflicts, or following conventional commit standards."
metadata:
  tags: git, version-control, branching, commits, collaboration
---

# Git Workflow

## Branch Management

**Naming conventions:**
- `feature/description` — new features
- `bugfix/description` — bug fixes
- `hotfix/description` — urgent fixes
- `refactor/description` — code refactoring
- `docs/description` — documentation updates

```bash
git checkout -b feature/feature-name
git push -u origin feature/feature-name
```

## Committing

**Conventional commit format:**
```
<type>(<scope>): <subject>

<body>

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
```

**Types:** `feat`, `fix`, `docs`, `style`, `refactor`, `test`, `chore`, `ci`

**sc_tools examples:**
- `feat(qc): add filter_spots() for visium_hd modality`
- `fix(ingest): correct load_imc_sample() processed_dir handling`
- `test(pp): add unit tests for run_scvi GPU detection`
- `ci: pin pulp<2.8 for snakemake 7 compatibility`

**Always use HEREDOC for multi-line messages:**
```bash
git commit -m "$(cat <<'EOF'
feat(qc): add generate_post_celltyping_report()

Implements the fourth QC report phase with biological metrics
(silhouette, cluster purity) meaningful only post-celltyping.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

## Staging

```bash
git add specific_file.py other_file.py   # prefer specific files over git add .
git diff --staged                         # review before committing
git add -p                               # interactive staging for partial changes
```

Never commit: `.env`, credentials, files >1MB, `*.h5ad`, `*.ai`, `*.tiff`.

## Merging and Conflict Resolution

```bash
git fetch origin
git merge origin/main                    # merge main into feature branch
# resolve conflicts, then:
git add resolved_file.py
git commit
```

For complex rebases: resolve conflicts one commit at a time, not all at once.

## sc_tools Safety Rules

- NEVER force push to main
- NEVER skip hooks (`--no-verify`) without explicit user request
- NEVER amend a published commit — create a new one
- Always check CI (`gh run list -L 1`) after pushing
- Stage specific files; avoid `git add .` to prevent committing large data files accidentally
