---
name: hpc-monitor-job
description: "[Tier 2 — Workflow] Submit SLURM jobs and monitor them to completion without blocking the conversation. Use when submitting sbatch jobs on HPC (brb, cayuga), running SpaceRanger, integration benchmarks, or any long-running cluster work."
tier: 2
---

# HPC Monitor Job

**Core principle:** Job submission and monitoring are one atomic workflow. Never submit a SLURM job and walk away — always attach a background monitor, verify outputs on completion, and update project docs.

**Announce at start:** "I'm using the hpc-monitor-job skill to submit and track this SLURM job."

## When to Use

**Use when:**
- Submitting any `sbatch` job on brb or cayuga
- Running SpaceRanger, scVI, deconvolution, or any long-running HPC task
- User asks to check on a previously submitted job
- Any pipeline phase that runs on SLURM

**Do NOT use when:**
- Running quick SSH commands (< 1 minute)
- Local-only work (no HPC involved)
- Dry-run or validation-only tasks

## The Pattern

### 1. Submit the Job

Write the sbatch script, copy to HPC, submit:
```bash
scp script.sbatch brb:/path/to/workdir/
ssh brb "cd /path/to/workdir && mkdir -p logs && sbatch script.sbatch"
```

Confirm submission with `squeue -j <job_id>`.

### 2. Start Background Monitor (MANDATORY)

Immediately launch a background Bash task using `run_in_background`:

```bash
# Template — adapt POLL_INTERVAL, OUTPUT_DIR, verification logic
while true; do
    STATUS=$(ssh -o ConnectTimeout=10 brb \
      "sacct -j ${JOB_ID} --format=JobID%-20,State%-12,ExitCode%-8,Elapsed%-12 --noheader 2>/dev/null" \
      2>/dev/null)

    # Print status
    echo "[$(date +%H:%M:%S)] $STATUS"

    # Check if all tasks finished
    PENDING=$(echo "$STATUS" | grep -cE 'RUNNING|PENDING' || true)
    if [[ "$PENDING" -eq 0 ]]; then
        # Verify outputs on HPC
        echo "[VERIFY] Checking expected outputs..."
        # ... verification logic ...
        exit 0
    fi

    sleep 90
done
```

**Key parameters:**
- Poll interval: **90 seconds** (Lustre is slow; faster polling wastes SSH connections)
- SSH timeout: `ConnectTimeout=10` to avoid hangs
- Use `sacct` (not `squeue`) — works after job finishes

### 3. Continue Working

Tell the user the background task ID. Continue with other work in the conversation. Do NOT block.

### 4. Check Results (when prompted or when background task completes)

Read background task output via `TaskOutput`. Then:

**If COMPLETED:**
- SSH to verify expected output files exist (binned_outputs, segmented_outputs, .h5ad, etc.)
- Update project Mission.md — mark tasks `[x]`, update status line
- Add Journal.md entry with date, job ID, results
- Update Journal.md with dated entry

**If FAILED:**
- Check SLURM .err log: `ssh brb "cat /path/to/logs/job_<id>.err | tail -50"`
- Check SpaceRanger/pipeline error files: `find <output_dir> -name "_errors" -exec cat {} \;`
- Report root cause to user
- Suggest fix or resubmission

### 5. Manifest and Doc Hygiene

After any Phase 0 job completes, also:
- Verify slide IDs match between manifest and CytAssist TIF filenames
- Update `metadata/phase0/all_samples.tsv` if corrections needed
- Update overall Phase 0a status count in Mission.md

## SSH and Lustre Rules

- **Always** use `ssh -o ConnectTimeout=10` in monitors (Lustre can hang indefinitely)
- **Prefer** `ls -d` or `find -maxdepth 1` over `ls` on large Lustre directories
- **Never** `cat` large files over SSH — copy to `/tmp` on the compute node first
- **Symlink FASTQs** rather than copy when source is on the same filesystem

## Array Job Pattern

For multiple samples, use SLURM arrays:
```bash
#SBATCH --array=0-N
SAMPLES=("sample1" "sample2" ...)
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
```

Monitor checks all array tasks before declaring completion.

## Integration with Other Skills

- **dispatching-parallel-agents:** Use one agent to submit, this skill to monitor. Multiple monitors can run concurrently for independent jobs.
- **journal-and-mission:** Always update docs on completion (this skill enforces it).
- **verification-before-completion:** Output verification is built into the monitor — no separate verification step needed for the HPC portion.

## Common Mistakes

- **Submitting and forgetting:** Always start the background monitor. This is the whole point.
- **Polling too fast:** 90s minimum. Faster wastes SSH connections and Lustre hates it.
- **Blocking the conversation:** Use `run_in_background`. Never `sleep` in the foreground.
- **Not verifying outputs:** A job can exit 0 but produce incomplete outputs (e.g. binned but no cell seg).
- **Not updating docs:** If you verified a job completed, update Mission.md in the same turn.
