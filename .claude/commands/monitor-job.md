Monitor a SLURM job on HPC until completion, then verify outputs.

## Inputs
- **$ARGUMENTS**: Either a SLURM job ID, or "last" to check the most recent background monitor task.

## Behavior

### If given a job ID:
1. SSH to brb and run `sacct -j <job_id> --format=JobID%-20,JobName%-30,State%-12,ExitCode%-8,Elapsed%-12,MaxRSS%-12 --noheader`
2. If any tasks are still RUNNING or PENDING, start a background monitor:
   - Use `run_in_background` Bash task
   - Poll `sacct` every 90 seconds
   - On completion: verify expected output files exist, print summary
   - Report failures with error details from SLURM logs
3. Tell the user the background task ID so they can check later
4. Continue the conversation — do NOT block waiting

### If "last" or no argument:
1. Check the most recent background monitor task via `TaskOutput` (non-blocking)
2. Report current status

### On job completion (when checking results):
1. Verify outputs exist on brb (check for expected files like `outs/`, `binned_outputs/`, `segmented_outputs/`)
2. If COMPLETED: update the relevant project Mission.md (mark tasks `[x]`) and Journal.md
3. If FAILED: SSH to check error logs (`_errors` files, SLURM .err logs), report root cause, suggest fix
4. Update any batch manifests or status files as needed

## Key rules
- **Never block the conversation** waiting for a job. Always use `run_in_background`.
- Poll interval: 90 seconds (not faster — Lustre is slow)
- Use `ssh -o ConnectTimeout=10 brb` to avoid SSH hangs
- After verification, always update Mission.md and Journal.md with results
- See skills.md section 18.7 for multi-agent orchestration patterns
