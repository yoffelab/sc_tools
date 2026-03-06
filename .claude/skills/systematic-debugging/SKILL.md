---
name: systematic-debugging
tier: 1
description: "[Tier 1 — Core] Mandatory four-phase debugging process: root cause first, no symptom fixes. Use whenever a bug, test failure, pipeline error, or unexpected behavior needs investigation."
---

# Systematic Debugging

**Core principle:** ALWAYS find root cause before attempting fixes. Symptom fixes are failure.

NO FIXES WITHOUT ROOT CAUSE INVESTIGATION FIRST.

## When to Use

Any time something is broken: test failures, SLURM job errors, pipeline crashes, unexpected outputs, import errors, OOM, wrong results.

## The Four Mandatory Phases

### Phase 1: Root Cause Investigation
- Read error messages fully — do not skim
- Reproduce the issue consistently before investigating
- Examine recent changes (git log, git diff)
- Gather diagnostic evidence across system components
- Trace data flow **backward** through call stacks

### Phase 2: Pattern Analysis
- Locate working examples of the same operation
- Compare against reference implementations completely
- Identify **all** differences between working and broken code
- Understand dependencies

### Phase 3: Hypothesis and Testing
- Form specific, testable hypotheses
- Make **one** minimal change at a time
- Verify each result before moving to the next
- Do not assume — if there is a gap in understanding, investigate it

### Phase 4: Implementation
- Write a failing test that reproduces the bug before fixing
- Implement a **single** targeted fix addressing the root cause
- Verify no regressions
- **If 3+ fix attempts fail: STOP and question the architecture** — do not continue symptom patching

## Critical Escalation Rule

Three failed fixes = structural problem. Do not continue. Redesign or ask for help.

## Red Flags — Restart the Process

- "Quick fix for now, investigate later"
- Proposing a solution before tracing data flow
- Attempting another fix after two prior failures
- "It works on my machine"
- Checking only recent code; ignoring dependencies

## sc_tools-specific notes

- For SLURM failures: check `sacct -j JOBID --format=State,ExitCode,MaxRSS` before re-submitting
- For h5ad errors: check file integrity with `h5py.File(path, 'r')` before loading via anndata
- For GPU OOM: check `nvidia-smi` and reduce `batch_size` before increasing `--mem`
- For conda/pip install hangs: switch to pip or mamba immediately (see skills.md §18.6)
