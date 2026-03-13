---
status: active
created: 2026-03-12
summary: "The repo has rules (CLAUDE"
---
wik# Plan: Claude Code Hooks + Project-Scoped Container Mounts

## Context

The repo has rules (CLAUDE.md, skills.md) about running code inside containers via `./scripts/run_container.sh` and using the `sc_tools` conda env, but nothing enforces these on the AI agent. The only existing hook is a PostToolUse lint check on `sc_tools/*.py` edits.

Additionally, `run_container.sh` currently mounts the entire repo as `/workspace` for all projects — meaning any project can read/write any other project. For modularity and safety, project containers should be scoped: project sees only its own dir + sc_tools (read-only).

**One conda env is sufficient** — projects differ in data and scripts, not dependencies. Platform-specific deps use optional extras (`pip install sc-tools[imc]`). Modularity comes from volume scoping, not separate environments.

**Goal:**
1. Add Claude Code hooks enforcing container usage and conda env
2. Modify `run_container.sh` to scope volume mounts per project vs sc_tools development

---

## Changes

### 1. `.claude/settings.json` — Add PreToolUse hook

Add a PreToolUse hook on Bash that:
- **Blocks** bare `python projects/...` or `python scripts/...` commands (exit 2) with a message suggesting the container equivalent
- **Warns** if `CONDA_DEFAULT_ENV != sc_tools` for host-side dev commands (`pytest`, `make`, `ruff`)
- **Allows** `python -c`, `python -m`, `make lint`, `pytest`, `ruff`, container commands

Hook script (inline):
```bash
cmd=$(echo "$CLAUDE_TOOL_INPUT" | python3 -c "import sys,json; print(json.load(sys.stdin).get('command',''))" 2>/dev/null)

# Block bare python on project/shared scripts
if echo "$cmd" | grep -qE '^python3?\s+projects/'; then
  proj=$(echo "$cmd" | grep -oE 'projects/[^/]+/[^/]+')
  echo "[hook] BLOCKED: Run project scripts inside a container:"
  echo "  ./scripts/run_container.sh $proj python <script>"
  exit 2
fi
if echo "$cmd" | grep -qE '^python3?\s+scripts/[^o]'; then
  echo "[hook] BLOCKED: Run shared scripts inside a container:"
  echo "  ./scripts/run_container.sh <project_path> python scripts/<script>"
  exit 2
fi

# Conda env warning for host-side dev commands
if echo "$cmd" | grep -qE '^(python|pytest|make)\b'; then
  if [ "$CONDA_DEFAULT_ENV" != "sc_tools" ] && [ -z "$APPTAINER_CONTAINER" ]; then
    echo "[hook] WARNING: conda env 'sc_tools' is not active. Run: conda activate sc_tools"
  fi
fi
```

### 2. `scripts/run_container.sh` — Project-scoped volume mounts

Modify the Docker and Apptainer runners to use different mount strategies based on the project path:

**sc_tools development** (project path is `.` or `sc_tools`):
- Mount full repo: `-v "$REPO_ROOT:/workspace"` (current behavior)
- Working directory: `/workspace`

**Project execution** (project path matches `projects/<modality>/<project>`):
- Mount project dir read-write: `-v "$PROJECT_PATH:/workspace/project"`
- Mount sc_tools read-only: `-v "$REPO_ROOT/sc_tools:/workspace/sc_tools:ro"`
- Mount scripts read-only: `-v "$REPO_ROOT/scripts:/workspace/scripts:ro"`
- Mount pyproject.toml etc read-only for imports to work
- Working directory: `/workspace/project`
- `PROJECT_DIR=/workspace/project`

Implementation in `run_docker()`:
```bash
run_docker() {
  local IMAGE="${SC_TOOLS_IMAGE:-sc_tools:latest}"
  # ... existing docker checks ...

  if [[ "$PROJECT_REL" == "." || "$PROJECT_REL" == sc_tools* ]]; then
    # sc_tools development: full repo access
    local OPTS=(-v "$REPO_ROOT:/workspace" -w "/workspace"
                -e "PROJECT_DIR=/workspace")
  else
    # Project execution: scoped mounts
    local OPTS=(
      -v "$PROJECT_PATH:/workspace/project"
      -v "$REPO_ROOT/sc_tools:/workspace/sc_tools:ro"
      -v "$REPO_ROOT/scripts:/workspace/scripts:ro"
      -v "$REPO_ROOT/pyproject.toml:/workspace/pyproject.toml:ro"
      -v "$REPO_ROOT/skills.md:/workspace/skills.md:ro"
      -w "/workspace/project"
      -e "PROJECT_DIR=/workspace/project"
    )
  fi

  if [[ ${#USER_CMD[@]} -eq 0 ]]; then
    exec docker run -it --rm "${OPTS[@]}" "$IMAGE" bash
  else
    exec docker run --rm "${OPTS[@]}" "$IMAGE" "${USER_CMD[@]}"
  fi
}
```

Same pattern for `run_apptainer()` using `--bind` with `:ro` suffix.

For `run_direct()` (no container): keep as-is (conda env handles isolation; no volume scoping possible without containers).

---

## Files to Modify

| File | Action | Description |
|------|--------|-------------|
| `.claude/settings.json` | **Modify** | Add PreToolUse Bash hook (container + conda enforcement); keep existing PostToolUse lint hook |
| `scripts/run_container.sh` | **Modify** | Add project-scoped volume mounts (project RW + sc_tools RO) vs full-repo mount for sc_tools dev |

---

## Verification

```bash
# 1. Hook blocks bare python on project scripts
python projects/visium/ggo_visium/scripts/score_gene_signatures.py
# → should be blocked with container suggestion

# 2. Hook blocks bare python on shared scripts
python scripts/run_qc_report.py
# → should be blocked

# 3. Hook allows dev commands
make lint        # → allowed
pytest sc_tools/tests/test_pp.py -v  # → allowed (warns if wrong conda env)
python -c "import sc_tools; print('ok')"  # → allowed

# 4. Container command allowed
./scripts/run_container.sh projects/visium/ggo_visium python scripts/score_gene_signatures.py
# → allowed, runs in container

# 5. Scoped mounts work for project
./scripts/run_container.sh projects/visium/ggo_visium bash -c "ls /workspace/"
# → should see project/ sc_tools/ scripts/ but NOT other projects

# 6. Full mount works for sc_tools dev
./scripts/run_container.sh . bash -c "ls /workspace/"
# → should see full repo

# 7. Existing lint hook still works
# Edit sc_tools/pp/normalize.py → should auto-lint
```
