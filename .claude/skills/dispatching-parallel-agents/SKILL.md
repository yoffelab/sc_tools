---
name: dispatching-parallel-agents
tier: 1
description: "[Tier 1 — Core] Handle 2+ independent problems simultaneously instead of sequentially. Use when facing multiple unrelated failures, independent pipeline phases, or embarrassingly parallel tasks (e.g. per-sample SLURM jobs, parallel integration methods)."
---

# Dispatching Parallel Agents

**Core principle:** Dispatch one agent per independent problem domain. Let them work concurrently.

## When to Use

**Use when:**
- 3+ test files failing with different root causes
- Multiple pipeline phases that are independent (`demographics` + `scoring`)
- Per-sample processing (SpaceRanger, QC, deconvolution per library)
- Each problem is self-contained with no shared state

**Do NOT use when:**
- Failures are related (fixing one might fix others)
- Agents would edit the same files
- Full system context is needed across all problems
- Shared state (same SLURM resources, same output dir)

## The Pattern

### 1. Identify Independent Domains
Group problems by what is broken or what needs to be done:
- File A tests: one root cause
- File B tests: different root cause
- Sample 1 processing: independent of Sample 2

### 2. Create Focused Agent Tasks
Each agent gets:
- **Specific scope:** one file, one sample, one method
- **Clear goal:** what to produce or fix
- **Constraints:** what NOT to touch
- **Expected output:** what to return

### 3. Dispatch in Parallel (single message, multiple Task tool calls)

```python
# Integration benchmark — 9 methods in parallel
for method in ["scvi", "harmony", "cytovi", "scvi+harmony", "raw", "pca", "bbknn", "scanorama", "desc"]:
    Task(agent="general-purpose",
         prompt=f"Run integration benchmark for method={method}. "
                f"Read results/adata.normalized.p3.h5ad. "
                f"Save to results/tmp/integration_test/{method}.h5ad.")
```

### 4. Review and Integrate
- Read each agent's summary
- Verify no file conflicts
- Run full test suite to confirm all changes work together

## Agent Prompt Structure

Good prompts are:
1. **Focused** — one clear problem domain
2. **Self-contained** — all context the agent needs
3. **Specific about output** — what to return

```markdown
Fix the 2 failing tests in sc_tools/tests/test_qc.py related to filter_spots():

1. test_filter_spots_min_counts — expects spots with total_counts < 500 removed
2. test_filter_spots_modality_visium — expects different defaults for visium

Do NOT change other functions. Return: what you found and what you fixed.
```

## sc_tools-specific parallel patterns

| Situation | Parallel split |
|-----------|---------------|
| SpaceRanger batch1 + batch2 | Two submission agents, one per batch |
| 9-method integration benchmark | 9 agents, one per method |
| `demographics` + `scoring` | Two analysis agents from same `preprocess` checkpoint |
| Deconvolution per library | One agent per library (GPU-limited: max 2 concurrent) |
| Per-sample QC (16 samples) | SLURM array (not agents); agents submit array jobs |

## Common Mistakes

**Too broad:** "Fix all the tests" — agent gets lost
**No constraints:** Agent may refactor unrelated code
**Vague output:** "Fix it" — you don't know what changed
**Shared state:** Two agents editing the same AnnData checkpoint
