---
name: verification-before-completion
tier: 1
description: "[Tier 1 — Core] Run fresh verification before claiming any task is done. Use before every success claim, commit, PR creation, or task transition. No 'should work', 'probably fixed', or 'Done!' without evidence."
---

# Verification Before Completion

**Core principle:** NO COMPLETION CLAIMS WITHOUT FRESH VERIFICATION EVIDENCE.

## The 5-Step Gate (mandatory before any success claim)

1. Identify which command validates the claim
2. Execute the complete command freshly — no cached results
3. Read the **full** output including exit codes
4. Verify output actually confirms the claim
5. Only then communicate the result with evidence

## Prohibited Patterns

Never say:
- "should work now"
- "probably fixed"
- "seems to be passing"
- "Great! Done!" before running verification
- "I've implemented X" without showing test/lint output

## What Counts as Verification

| Claim | Required evidence |
|-------|------------------|
| "Tests pass" | Fresh `pytest` run with zero failures shown |
| "Lint clean" | `ruff check sc_tools` with no output |
| "CI passes" | `gh run list -L 1` showing `completed success` |
| "Script runs" | Actual output of running the script |
| "File saved correctly" | Read the file back and confirm contents |
| "SLURM job finished" | `sacct -j JOBID` showing `COMPLETED` state |

## sc_tools-specific checkpoints

After every commit:
```bash
gh run list -L 1   # wait for CI, check completed success
```

After every code change:
```bash
make lint          # ruff check + ruff format --check
pytest sc_tools/tests/ -v --tb=short -q
```

After SLURM job submission:
```bash
squeue -u juk4007  # confirm job is queued/running
sacct -j JOBID --format=JobID,State,ExitCode,Elapsed,MaxRSS  # after completion
```

## Key Insight

Past unverified claims have led to undefined functions shipping and incomplete features reaching production. Fresh evidence is non-negotiable — previous runs, code reading, or agent reports are insufficient.
