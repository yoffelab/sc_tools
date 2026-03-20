# Pitfalls Research

**Domain:** Agent-native CLI for computational biology pipelines (single-cell/spatial)
**Researched:** 2026-03-20
**Confidence:** HIGH (grounded in codebase evidence + domain research)

## Critical Pitfalls

### Pitfall 1: CLI Output as Unstable Contract

**What goes wrong:**
The CLI returns JSON to agents, but field names, nesting, or types change between versions. Agents that parsed `{"cells": 50000}` break silently when it becomes `{"n_obs": 50000}`. Since agents adopt features atomically (not gradually like humans), a single field rename breaks every agent workflow simultaneously.

**Why it happens:**
Developers treat CLI output as informal logging rather than a versioned API. During early development, fields get renamed for "consistency" without realizing downstream consumers have already hardcoded the old names. sc_tools already has this pattern -- dual phase nomenclature (p1 vs qc_filter) coexists, and the registry schema has changed 7 times in recent migrations.

**How to avoid:**
- Define a JSON output schema (TypedDict or Pydantic model) for each command from day one
- Never remove or rename a field -- add new fields, deprecate old ones with a version flag
- CI test that validates CLI output against the schema on every commit
- Version the output format: `{"version": "1.0", "data": {...}}`

**Warning signs:**
- Multiple output format changes in the same sprint
- Agent code that uses `result.get("field", result.get("old_field"))` fallback patterns
- No schema file or TypedDict for CLI output structures

**Phase to address:**
Phase 1 (CLI foundation) -- schema must be defined before any command ships

---

### Pitfall 2: Schema Migration Ratchet (Already Happening)

**What goes wrong:**
The registry database schema keeps changing because the domain model is not yet stable. Each migration adds complexity: old tables coexist with new ones, dual-write patterns emerge, rollbacks become untested, and eventually the migration chain itself becomes a blocker. sc_tools has 7 unapplied migrations (0012-0018), a four-layer schema that creates 7 new tables with cross-join data migration, and a dual-write pattern in `register_dataset()` that writes to both old and new tables.

**Why it happens:**
Schema changes feel cheap in SQLAlchemy/Alembic -- `alembic revision --autogenerate` creates a file, and the developer moves on. But each migration is a permanent node in a directed graph. After 18 migrations with some creating tables that later migrations drop, the upgrade-from-scratch path becomes fragile. The four-layer migration (0013) contains cross-join logic that could timeout on large registries.

**How to avoid:**
- **Do not use the database for provenance in Phase 1.** The PROJECT.md already calls for file-based provenance (JSON sidecars) first. Honor this decision strictly.
- Freeze the existing DB schema -- no new migrations until file-based provenance proves what fields are actually needed
- When DB provenance eventually lands, squash all migrations into a single baseline before adding new ones
- Never use ORM models inside migration scripts (they drift from the migration's intent)

**Warning signs:**
- More than 2 migration files created in a single sprint
- Any migration that both creates AND drops tables
- Dual-write patterns in application code
- Migration tests that are skipped (sc_tools already has this: `test_migrations.py` line 118 is a `pass`)

**Phase to address:**
Phase 1 must use file-based provenance exclusively. DB migration cleanup is a separate, later phase after the domain model stabilizes through file-based usage.

---

### Pitfall 3: Loading Full AnnData When CLI Only Needs Metadata

**What goes wrong:**
A CLI command like `sct info sample.h5ad` loads the entire 25G file into memory just to report obs count and var names. AnnData in memory uses roughly 4x disk size, so a 25G file needs ~100G RAM. The CLI becomes unusable on login nodes or laptops for basic queries. This already happened: the PROJECT.md origin story describes an agent-written benchmark script that OOM'd at 44G.

**Why it happens:**
`scanpy.read_h5ad()` loads everything by default. Developers use `adata = sc.read_h5ad(path)` out of habit. Even backed mode (`backed='r'`) still loads obs, var, uns, and obsp into memory -- only X stays on disk. For a 2.5M-cell dataset, obs alone can be several GB.

**How to avoid:**
- Use `h5py` directly for metadata-only commands: `h5py.File(path)['obs'].attrs` gives obs columns without loading data
- Categorize CLI commands into tiers: metadata-only (h5py), summary (backed mode + subsample), compute (full load with memory budget)
- Set explicit memory budgets per command tier and fail fast with a clear message if the file exceeds the budget
- For benchmark reports, load embeddings only (not raw counts) via h5py dataset slicing

**Warning signs:**
- CLI commands that import scanpy at the top level (forces scipy/numpy overhead even for help text)
- Any command that calls `sc.read_h5ad()` without `backed=True` on files over 1G
- CLI startup time exceeding 2 seconds (indicates heavy imports)

**Phase to address:**
Phase 1 (CLI foundation) -- command dispatch must use lazy imports; Phase 2 (data commands) must implement the tiered loading strategy

---

### Pitfall 4: Interactive Prompts and Unbounded Output Kill Agent Workflows

**What goes wrong:**
A CLI command pauses for confirmation ("Delete 50 samples? [y/N]"), pages output through `less`, or prints a progress bar to stderr. The agent's workflow hangs indefinitely on the prompt, the pager, or drowns in spinner characters. AWS CLI v2 broke thousands of CI jobs in 2019 by adding an interactive pager as default behavior.

**Why it happens:**
Human-friendly defaults (confirmations, pagers, color codes) are hostile to programmatic consumers. Developers add a confirmation prompt for destructive operations and forget that agents cannot type "y". Progress bars seem harmless but pollute stderr parsing.

**How to avoid:**
- Never prompt for input by default. Destructive operations require `--confirm` flag (not interactive prompt)
- Detect `NO_COLOR` and `CI` environment variables; suppress color/spinners when set
- All output goes to stdout (structured JSON). Diagnostics go to stderr. Never mix them
- Cap JSON output at a sensible default (e.g., 100 items) with `--limit` and pagination cursors
- Test every command in headless mode (pipe stdout to /dev/null, ensure no hangs)

**Warning signs:**
- Any `input()` or `click.confirm()` call in CLI code
- Output that changes based on terminal width
- Progress bars written to stdout instead of stderr
- Commands that produce >10KB of output without pagination

**Phase to address:**
Phase 1 (CLI foundation) -- non-interactive-by-default must be an architectural invariant, not a feature added later

---

### Pitfall 5: Observation ID Collisions in Multi-Modal Assembly

**What goes wrong:**
When assembling MuData from independently processed modalities (Visium + scRNA-seq for the same patient), observation IDs collide or silently fail to link. MuData treats observations with identical names across modalities as the same observation. If Visium uses barcodes (`ACGT-1`) and scRNA-seq uses different barcodes for the same patient, the modalities appear completely disjoint. Conversely, if both use sequential integers (`0, 1, 2...`), unrelated cells get falsely linked.

**Why it happens:**
Each modality is processed independently with its own ID scheme. Spatial platforms use spatial barcodes, scRNA-seq uses library-specific barcodes, IMC uses ROI-cell combinations. There is no natural shared identifier at the cell level -- linkage happens at the patient/sample level, not cell level. Developers assume MuData handles this automatically, but it only handles the container structure, not the semantic linking.

**How to avoid:**
- Define a patient/sample ID contract early: `{project}_{patient}_{sample}_{modality}` as the canonical obs index prefix
- Assembly commands must validate that modality obs indices do NOT accidentally intersect (unless truly paired measurements)
- Store linkage metadata in MuData.obs (patient_id, sample_id) rather than relying on index matching
- Implement a `sct mudata validate` command that checks for accidental index collisions and missing linkage fields

**Warning signs:**
- MuData.obs has fewer rows than expected (silent intersection)
- Integration methods fail with "no shared observations" errors
- Patient-level queries return cells from only one modality

**Phase to address:**
Late phase (MuData assembly) -- but the patient/sample ID contract must be defined in Phase 1 metadata design

---

### Pitfall 6: Lazy Import Tax vs. Startup Speed

**What goes wrong:**
Every CLI command pays a 3-8 second startup tax because `import sc_tools` triggers `import scanpy`, which imports numpy, scipy, matplotlib, and anndata. Even `sct --help` takes 5 seconds. Agents calling the CLI 50 times per workflow spend 4 minutes just on imports. This makes the CLI unusable for the rapid chained-command patterns agents prefer.

**Why it happens:**
Python's import system is eager. `from sc_tools.pp import normalize` triggers the entire sc_tools package init. Scientific Python stacks are particularly heavy -- scanpy alone imports ~30 packages.

**How to avoid:**
- CLI entrypoint must be a thin dispatcher that imports only argparse/sys/json
- Each subcommand module uses function-level imports: `def run_qc(): import scanpy as sc`
- Never import scanpy, anndata, or matplotlib at module level in CLI code
- Measure startup time in CI: `time sct --help` must complete in <0.5s
- Consider a long-running daemon mode (`sct daemon`) for agent workflows that amortizes import cost

**Warning signs:**
- `sct --help` takes >1 second
- Adding a new CLI command increases import time for all commands
- `__init__.py` files that import submodules eagerly

**Phase to address:**
Phase 1 (CLI foundation) -- the entrypoint architecture must enforce lazy imports from the start

---

## Technical Debt Patterns

Shortcuts that seem reasonable but create long-term problems.

| Shortcut | Immediate Benefit | Long-term Cost | When Acceptable |
|----------|-------------------|----------------|-----------------|
| Wrapping existing scripts as CLI commands verbatim | Ship fast | Inconsistent argument names, no structured output, different error patterns per command | Never -- refactor the interface even if internals stay the same |
| Storing provenance in AnnData .uns | No external files needed | uns is untyped dict, no schema validation, lost on subset/concat, bloats file size | Only for ephemeral run metadata (timestamp, sc_tools version) |
| SQLite for concurrent agent access | Zero-config setup | Write locks cause timeout under multi-agent load; already identified in CONCERNS.md | Only for single-user local dev; document Postgres for production |
| Printing Python tracebacks as error output | Easy debugging | Agents parse tracebacks as "unknown error", wastes tokens, exposes internals | Never in CLI mode; tracebacks go to --debug log file only |
| `click.Group` with implicit command discovery | Auto-registers new commands | Import-time side effects, hard to control loading order, all commands imported on startup | Never -- use explicit command registration with lazy loading |

## Integration Gotchas

Common mistakes when connecting CLI to existing sc_tools internals.

| Integration | Common Mistake | Correct Approach |
|-------------|----------------|------------------|
| sc_tools.qc.report | Passing AnnData objects through CLI (serialization overhead) | Pass file paths; let the report function load what it needs |
| sc_tools.pp.recipes | Assuming GPU availability | CLI must detect backend (`has_rapids()`) and report it in structured output; never silently fall back without telling the agent |
| sc_tools.registry | Calling registry methods that assume SQLAlchemy session exists | CLI provenance layer must be file-based (JSON sidecars), not registry-dependent |
| sc_tools.mcp | Duplicating MCP tool logic in CLI commands | Share a common implementation; CLI and MCP are both thin interfaces over the same function |
| Snakemake rules | CLI commands that duplicate Snakemake rule logic | CLI commands should be callable FROM Snakemake rules, not compete with them |
| h5ad file locking | Multiple CLI commands writing to the same h5ad | h5ad is not concurrent-write-safe; use output-to-new-file pattern, never in-place modification |

## Performance Traps

Patterns that work at small scale but fail as data grows.

| Trap | Symptoms | Prevention | When It Breaks |
|------|----------|------------|----------------|
| Loading full AnnData for summary stats | OOM on login nodes, 30s+ command time | Use h5py for metadata, backed mode for summaries, full load only for compute | >5G files (~500K cells) |
| Computing all benchmark metrics in one pass | Memory spike from holding multiple integration results | Stream metrics: load one method at a time, compute, write, release | >3 integration methods on >1M cells |
| JSON sidecar files with embedded arrays | Slow parse, large files, git diff noise | Store arrays as references to h5ad datasets, not inline JSON | >10K observations in provenance |
| Unbounded `sct list` output | Agent context overflow, slow parsing | Default limit of 50 items, cursor-based pagination | >100 registered datasets |
| Eager validation of all inputs before any work | Long delay before first output, no partial results | Validate per-step, emit structured progress as each step completes | Pipelines with >5 steps |
| k-NN on full dataset for QC | OOM; current code already subsamples to 50K (FIT_LIMIT) | Document the subsampling; expose as CLI parameter | >50K spatial spots |

## Security Mistakes

Domain-specific security issues.

| Mistake | Risk | Prevention |
|---------|------|------------|
| DB URL with credentials in SC_TOOLS_REGISTRY_URL | Credential leak in error messages, logs, agent output | Mask credentials in all output; use separate auth env vars; validate URL does not contain passwords |
| Clinical data (patient IDs, diagnoses) in CLI output | PHI exposure in agent logs or version control | Never include patient identifiers in default CLI output; require `--include-phi` flag with warning |
| h5ad files with embedded clinical metadata | PHI in shared analysis files | Validate that shared h5ad files have clinical columns stripped; add `sct scrub` command |
| Unvalidated file paths from agent input | Path traversal if CLI constructs paths from agent-provided strings | Resolve and validate all paths against project root; reject `../` patterns |

## UX Pitfalls

Mistakes that make the CLI frustrating for both agents and humans.

| Pitfall | User Impact | Better Approach |
|---------|-------------|-----------------|
| Inconsistent verb naming (run_qc vs generate_report vs do_integration) | Agents cannot predict command names | Consistent verb set: `ingest`, `qc`, `preprocess`, `integrate`, `annotate`, `benchmark`, `assemble` |
| Different flags for the same concept (`--input`, `--adata`, `--h5ad`, `--file`) | Agents must learn per-command flag names | Single convention: `--input` for primary input, `--output` for primary output, always |
| Error messages that say "failed" without saying why or what to fix | Agent retries the same command, wastes compute | Structured error: `{"error": "file_not_found", "path": "/x/y.h5ad", "suggestion": "check path exists"}` |
| Commands that succeed silently (exit 0, no stdout) | Agent cannot confirm success | Always emit result JSON: `{"status": "ok", "output_path": "/x/y.h5ad", "n_obs": 50000}` |
| `--verbose` as only output control | Either too little or too much | Three levels: default (JSON summary), `--verbose` (JSON detailed), `--debug` (human-readable diagnostics to stderr) |
| Help text that lists flags but not examples | Agent must guess correct invocation | Every command help includes at least one complete working example |

## "Looks Done But Isn't" Checklist

Things that appear complete but are missing critical pieces.

- [ ] **CLI command "works":** Often missing structured error output -- verify the command returns valid JSON for ALL failure modes, not just the happy path
- [ ] **Provenance tracking "records" runs:** Often missing the input file hash -- verify that re-running with identical inputs produces identical provenance (deterministic)
- [ ] **JSON output "has all fields":** Often missing units and types -- verify that numeric fields include units (`size_mb`, not `size`) and timestamps are ISO 8601
- [ ] **Multi-modal assembly "creates MuData":** Often missing validation that modalities actually linked -- verify `.obs` has patient_id column and modality counts match expectations
- [ ] **`sct help` "lists commands":** Often missing argument descriptions -- verify each command's help includes parameter types, defaults, and one example
- [ ] **Exit codes "indicate success/failure":** Often only 0/1 -- verify at least: 0=success, 1=user error (bad args), 2=data error (corrupt file), 3=runtime error (OOM), 4=dependency missing
- [ ] **CLI "handles large files":** Often tested only on toy data -- verify each data command on a real 2.5M-cell file; measure peak RSS with `/usr/bin/time -v`
- [ ] **Schema migration "works":** Often tested only forward -- verify downgrade path; verify upgrade from clean DB; verify upgrade from every historical state

## Recovery Strategies

When pitfalls occur despite prevention, how to recover.

| Pitfall | Recovery Cost | Recovery Steps |
|---------|---------------|----------------|
| Output schema break | MEDIUM | Add back old field as alias; emit deprecation warning; bump minor version; give agents 2 releases to migrate |
| Migration chain corruption | HIGH | Dump current DB data as JSON; drop and recreate from squashed baseline; reimport data; never attempt to fix the chain in-place |
| OOM in CLI command | LOW | Kill process; re-run with `--subsample N` flag; file a bug if default should have subsampled |
| Agent workflow hung on prompt | LOW | Kill process; add `--yes` flag; add CI test for headless mode |
| MuData obs ID collision | HIGH | Rebuild MuData with corrected obs indices; audit all downstream analyses that used the incorrect assembly; no shortcut |
| Provenance sidecar drift (sidecar does not match h5ad) | MEDIUM | Re-derive provenance from h5ad metadata (uns, obs columns); add `sct provenance verify` command |

## Pitfall-to-Phase Mapping

How roadmap phases should address these pitfalls.

| Pitfall | Prevention Phase | Verification |
|---------|------------------|--------------|
| Output schema instability | Phase 1: CLI foundation | CI test validates all command outputs against TypedDict schemas |
| Schema migration ratchet | Phase 1: File-based provenance (avoid DB) | Zero new Alembic migrations created during CLI development phases |
| Full AnnData loading for metadata | Phase 2: Data commands | `sct info` on 25G file completes in <5s with <500MB RSS |
| Interactive prompts/unbounded output | Phase 1: CLI foundation | CI test runs every command with stdin=/dev/null, stdout piped; no hangs |
| Obs ID collision in MuData | Phase N: Multi-modal assembly | `sct mudata validate` command checks for collision; integration test with 2+ modalities |
| Import tax / slow startup | Phase 1: CLI entrypoint | `time sct --help` < 0.5s in CI |
| Inconsistent command naming | Phase 1: CLI design spec | Lint rule or test that all commands use approved verb set |
| PHI in CLI output | Phase 1: Output contract | No patient_id or diagnosis fields in default output; grep test on all output schemas |
| Provenance over-engineering | Phase 1-2: JSON sidecars | Provenance is a single JSON file per output; no DB writes for provenance in early phases |
| Agent error parsing failure | Phase 1: Error contract | Every CLI error path returns `{"error": "<type>", "message": "<text>"}` JSON |

## Sources

- [Anthropic: Writing effective tools for AI agents](https://www.anthropic.com/engineering/writing-tools-for-agents) -- tool design principles for agent consumption
- [InfoQ: Keep the Terminal Relevant: Patterns for AI Agent Driven CLIs](https://www.infoq.com/articles/ai-agent-cli/) -- output contracts, idempotency, discovery patterns
- [DEV: Writing CLI Tools That AI Agents Actually Want to Use](https://dev.to/uenyioha/writing-cli-tools-that-ai-agents-actually-want-to-use-39no) -- exit codes, structured output, non-interactive design
- [scanpy Issue #2365: backed mode OOM on large datasets](https://github.com/scverse/scanpy/issues/2365) -- memory pitfalls with large h5ad files
- [scanpy Issue #434: backed='r' does not reduce memory](https://github.com/scverse/scanpy/issues/434) -- backed mode loads obs/var/uns into memory
- [MuData documentation: multimodal data objects](https://mudata.readthedocs.io/en/latest/io/mudata.html) -- obs identity semantics across modalities
- [PingCAP: Best Practices for Alembic Schema Migration](https://www.pingcap.com/article/best-practices-alembic-schema-migration/) -- migration chain management
- [Alembic Discussion #1259: managing large sets of migration files](https://github.com/sqlalchemy/alembic/discussions/1259) -- migration squashing, performance
- [Tracer: Bioinformatics Pipeline Frameworks 2025](https://www.tracer.cloud/resources/bioinformatics-pipeline-frameworks-2025) -- provenance over-engineering lessons
- sc_tools codebase: `.planning/codebase/CONCERNS.md` -- existing migration debt, registry concurrency, dual-write patterns
- sc_tools codebase: `.planning/PROJECT.md` -- project constraints, data scale, decision to use file-based provenance first

---
*Pitfalls research for: Agent-native CLI for computational biology pipelines*
*Researched: 2026-03-20*
