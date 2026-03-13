---
status: in_progress
created: 2026-03-12
summary: Rsync sc_tools (~2.5 TB) from local Mac to cayuga HPC at ~/elementolab/sc_tools/
---

# Plan: Rsync sc_tools to cayuga

## Goal

Transfer the full `sc_tools` repository (~2.5 TB) from
`/Users/junbumkim/Documents/sc_tools/` to `cayuga:~/elementolab/sc_tools/`
using a 5-pass size-tiered rsync with a priority pre-pass for the most active projects.

---

## Infrastructure

| Item | Detail |
|------|--------|
| Script | `~/sync_sc_tools.sh` |
| Log | `~/sync_sc_tools.log` |
| tmux session | `sc_tools_sync` |
| Remote dir | `cayuga:~/elementolab/sc_tools/` (chmod 700) |
| rsync binary | `/opt/homebrew/bin/rsync` (v3.4.1 — required for `--append-verify`) |

**Resume command:** `tmux send-keys -t sc_tools_sync "~/sync_sc_tools.sh <pass>" ENTER`

---

## Pass Plan

| Pass | Arg | Content | Size est. | Status |
|------|-----|---------|-----------|--------|
| priority — robin | `priority` | `projects/visium_hd/robin/` full | ~200 GB | **DONE** (2026-03-04 19:19) |
| priority — robin_round2 | `priority` | `projects/visium_hd/robin_round2/` full | ~200 GB | **Partial** — interrupted multiple times; resume from `priority` |
| priority — ggo_visium | `priority` | `projects/visium/ggo_visium/` full | ~100 GB | Not started |
| Pass 1 | `1` | Code, configs, metadata (no data/results/figures/outputs) | ~5 MB | Not started (earlier attempts failed — old rsync lacked `--append-verify`) |
| Pass 2 | `2` | `--max-size=100m`, all dirs | ~1 GB | Not started |
| Pass 3 | `3` | `--max-size=1g` | ~50 GB | Not started |
| Pass 4 | `4` | `--max-size=10g` | ~300 GB | Not started |
| Pass 5 | `5` | No size limit (fastq tars, BAMs, git objects) | ~2+ TB | Not started |

---

## Current State (2026-03-12)

- **robin**: DONE
- **robin_round2**: Partially transferred; was mid-transfer on data tarballs when last paused (Mar 7). Shell became unresponsive after Ctrl+C — session may need to be re-checked.
- **ggo_visium, passes 1–5**: Not started

### Known issues

- Shell/tmux became unresponsive after Ctrl+C on Mar 7 — verify `tmux ls` and `tail ~/sync_sc_tools.log` before resuming
- VPN (WCM/Cornell) must be active before any rsync or SSH to cayuga
- Large `.git` pack files in `projects/imc/aapc/` will be slow in Pass 5

---

## Next Steps

1. **Verify tmux is alive**: `tmux ls` — if dead, just run `tmux new-session -d -s sc_tools_sync`
2. **Check log**: `tail -10 ~/sync_sc_tools.log` — confirm last state
3. **Resume**: `tmux send-keys -t sc_tools_sync "~/sync_sc_tools.sh priority" ENTER`
   - Will skip robin (already done), resume robin_round2, then auto-advance through ggo_visium → passes 1–5
4. **Monitor**: `tail -f ~/sync_sc_tools.log` or `tmux attach -t sc_tools_sync` (detach: Ctrl+B D)
5. **Verify remote after all passes**:
   ```bash
   ssh cayuga "du -sh ~/elementolab/sc_tools/"           # should be ~2.5 TB
   ssh cayuga "ls ~/elementolab/sc_tools/projects/"      # check all projects present
   ```

---

## Pause / Resume Reference

| Action | Command |
|--------|---------|
| Pause | `tmux send-keys -t sc_tools_sync C-c` |
| Resume from priority | `~/sync_sc_tools.sh priority` |
| Resume from pass N | `~/sync_sc_tools.sh N` |
| Check progress | `tail -30 ~/sync_sc_tools.log` |
| Watch live | `tmux attach -t sc_tools_sync` (detach: Ctrl+B D) |
| Check remote size | `ssh cayuga "du -sh ~/elementolab/sc_tools/"` |
