---
name: gh-cli
tier: 2
description: "[Tier 2 — Workflow] GitHub CLI reference for CI monitoring, PRs, issues, Actions, and releases. Use when checking CI status after a push, creating PRs, monitoring workflow runs, or managing GitHub resources from the command line."
allowed-tools: Bash
---

# GitHub CLI (gh)

## Most-Used Commands for sc_tools

### CI Monitoring (check after every push)
```bash
gh run list -L 1                          # latest run status
gh run list -L 5                          # last 5 runs
gh run watch                              # stream live output of latest run
gh run view RUN_ID                        # details of specific run
gh run view RUN_ID --log                  # full logs
gh run view RUN_ID --log-failed           # only failed step logs

# Wait for CI and check result
gh run watch && gh run list -L 1
```

### Pull Requests
```bash
gh pr create --title "title" --body "$(cat <<'EOF'
## Summary
- bullet points

## Test plan
- [ ] CI passes
- [ ] lint clean

🤖 Generated with Claude Code
EOF
)"

gh pr list                    # list open PRs
gh pr view 123                # view PR details
gh pr checks 123              # check CI status on PR
gh pr merge 123 --squash      # merge PR
```

### Workflow Control
```bash
gh workflow list              # list all workflows
gh workflow run ci.yml        # trigger workflow_dispatch
gh run cancel RUN_ID          # cancel a running job
gh run rerun RUN_ID           # rerun failed run
gh run rerun RUN_ID --failed-only  # rerun only failed jobs
```

### Releases
```bash
gh release list
gh release create v1.2.3 --title "v1.2.3" --notes "Release notes"
gh release view v1.2.3
```

### Issues
```bash
gh issue create --title "Bug: ..." --body "..."
gh issue list
gh issue close 42
```

## sc_tools-specific CI workflow

After every `git push`:
1. `gh run list -L 1` — confirm run started
2. Wait or `gh run watch`
3. If failed: `gh run view RUN_ID --log-failed` to see which job
4. Fix, push again, repeat

Jobs to check: Lint, Test (py3.11/3.12/3.13/3.14 × ubuntu/macos), Docs, Snakemake dry-run, Docker (main only).

## Output Formatting
```bash
# JSON output with jq
gh run list --json status,conclusion,name -L 5 | jq '.[] | {name, status, conclusion}'

# Filter only failed runs
gh run list -L 20 --json conclusion,name | jq '.[] | select(.conclusion=="failure")'
```

## Authentication
```bash
gh auth status        # check current auth
gh auth login         # re-authenticate
```
