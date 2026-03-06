---
name: code-review
tier: 2
description: "[Tier 2 — Workflow] Conduct thorough code reviews covering correctness, security, performance, and testing. Use when reviewing PRs, auditing sc_tools modules, or checking implementation quality."
allowed-tools: Read Grep Glob
---

# Code Review

## When to Use
- Reviewing pull requests or diffs
- Auditing sc_tools module implementations
- Checking security and performance of analysis scripts
- Verifying test coverage before merging

## Review Process (8 Steps)

### 1. Understand context
What does this change do? What phase of the pipeline does it touch?

### 2. Architecture assessment
- Does it follow sc_tools conventions (checkpoint names, obsm keys, project paths)?
- Is logic in the right layer (sc_tools vs project script)?
- No project-specific code in sc_tools; no root-level results/figures/metadata

### 3. Code quality
- Functions focused and <50 lines where possible
- Clear naming (no abbreviations that aren't standard in the domain)
- No magic numbers — constants should be named
- No deep nesting (>3 levels = refactor signal)

### 4. Security
- No hardcoded credentials or paths
- No shell injection in subprocess calls
- Validate inputs at system boundaries (user input, external files)

### 5. Performance
- No redundant `adata.copy()` calls on large AnnDatas
- No loading full h5ad when only `obs` is needed
- Deconvolution batched by `library_id`
- Large data ops use backed AnnData where appropriate

### 6. Testing
- New functions have unit tests in `sc_tools/tests/`
- Tests use fixtures from `conftest.py`, not real data
- Edge cases covered: empty AnnData, missing obs columns, wrong modality

### 7. Documentation
- Public functions have docstrings
- Non-obvious logic has inline comments
- Checkpoint inputs/outputs documented

### 8. Feedback
Categorize issues:
- **Critical** — correctness bugs, data corruption risk, security — must fix before merge
- **Important** — performance, missing tests, convention violations — fix before merge
- **Nice-to-have** — style, minor clarity — optional

## sc_tools-specific Checklist

- [ ] Checkpoint files use standard names (Architecture.md §2.1)
- [ ] `obs` keys match required metadata contracts (Architecture.md §2.2)
- [ ] Signature scores in `obsm['signature_score']`, not `obs`
- [ ] Colors in `uns['{name}_colors']`, not overwritten if length matches
- [ ] No files >1MB added to git
- [ ] `make lint` passes (ruff check + format)
- [ ] Tests added or updated
- [ ] CI passes after push
