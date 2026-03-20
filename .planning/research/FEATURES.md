# Feature Landscape

**Domain:** Agent-native CLI for computational biology (single-cell/spatial analysis)
**Researched:** 2026-03-20

## Table Stakes

Features users (agents and humans) expect. Missing = agents write throwaway scripts instead of using the CLI.

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| **Structured JSON output on all commands** | Agents parse JSON, not human tables. gh CLI, kubectl, and docker all do this with `--json` or `-o json`. Without it, agents regex-parse stdout -- brittle and error-prone. | Low | Default to JSON; `--human` flag for readable output. NDJSON for streaming large results. |
| **Self-describing help (`sct describe`)** | Agents need to discover what commands exist and what params they accept at runtime, not from training data. CLI-Anything and gh CLI both expose this. | Low | `sct list-commands --json` returns command names, params, types, defaults. `sct describe <cmd>` returns JSON schema for a single command. |
| **Semantic exit codes** | Exit 0 = success, 1 = user error (bad args), 2 = data error (validation failed), 3+ = application errors. Agents use exit codes for branching logic. Without them, agents must parse error messages. | Low | Document exit code contract. Keep stable across versions. |
| **Structured error reporting** | Errors as JSON with `{"error": "...", "code": "...", "suggestion": "..."}` to stderr. Agents need machine-readable errors to decide retry vs. abort vs. fix. | Low | Always JSON to stderr. Include actionable suggestion field. |
| **Verb-noun command structure** | `sct qc run`, `sct preprocess run`, `sct integrate run` -- predictable grammar agents can reason about. Like `kubectl get`, `kubectl apply`. | Low | Verb-first: `run`, `list`, `describe`, `validate`, `report`. Noun second: `qc`, `preprocess`, `integrate`, `benchmark`, `celltype`. |
| **Checkpoint validation** | `sct validate <phase> <file>` -- agents need to verify outputs before proceeding to the next phase. Existing `validate.py` already does this. | Low | Wrap existing `validate_checkpoint()`. Return JSON with pass/fail + specific issues. |
| **Non-interactive execution** | No interactive prompts. All params via flags, env vars, or config files. Agents cannot respond to prompts. | Low | `--yes` for confirmations. Fail fast if required params missing (never prompt). |
| **Idempotent operations** | Re-running a command with same inputs produces same outputs. Agents retry on failure -- if the tool mutates state on partial failure, agents can't recover. | Medium | Design commands as pure functions: inputs -> outputs. Don't modify input files in-place. |
| **Dry-run mode** | `sct preprocess run --dry-run` validates inputs and reports what would happen without executing. Agents use this for planning and cost estimation (runtime, memory). | Medium | Return JSON with estimated runtime, memory, output paths. Critical for 2.5M-cell datasets. |
| **Phase-aware pipeline status** | `sct status <project>` returns which phases are complete, what's next, what files exist. Agents need state awareness to decide next action. Wraps existing registry/pipeline DAG. | Low | JSON output: completed phases, available next phases, checkpoint file paths, validation status. |

## Differentiators

Features that set this CLI apart from scanpy-scripts or ad-hoc Python wrappers. Not expected, but make agent workflows dramatically more reliable.

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| **JSON sidecar provenance** | Every output file gets a `.provenance.json` sidecar recording: command, params, input files (with checksums), sc_tools version, timestamp, runtime, peak memory. Nextflow 25.04 added native lineage -- this is the file-based equivalent for non-workflow-manager usage. | Medium | File-based first (PROJECT.md decision). Enables reproducibility without a database. Agents can read provenance to understand what produced a file. |
| **Memory estimation before execution** | `sct preprocess estimate <file>` returns estimated peak memory based on cell count, gene count, and method. Prevents OOM on 2.5M-cell Visium HD datasets -- the exact problem that motivated this project. | Medium | Use heuristics: scVI ~ 4x dataset size, Harmony ~ 2x, PCA ~ 1.5x. Return JSON with recommendation (GPU vs CPU, subsampling threshold). |
| **Agent-oriented error taxonomy** | Errors categorized as `retryable` (OOM, timeout), `fixable` (bad param, missing dep), `fatal` (corrupt data). Agents use this to decide retry, fix, or escalate. Goes beyond exit codes. | Low | JSON error objects: `{"category": "retryable", "code": "OOM", "suggestion": "Add --subsample 500000 or use GPU"}`. |
| **Pre-computed benchmark ingestion** | `sct benchmark load --method harmony --file harmony.h5ad --method scvi --file scvi.h5ad` loads separate per-method h5ad files via h5py (not loading full datasets into memory), computes metrics, generates comparison report. Solves the 44G OOM problem. | High | h5py-backed selective reads. Subsampling for metrics. Core value proposition from PROJECT.md. |
| **MuData assembly from independent modalities** | `sct assemble --rna adata_rna.h5ad --spatial adata_visium.h5ad --patient-key patient_id` creates MuData linking modalities by patient/subject. No other comp bio CLI does this as a single command. | High | Uses muon/MuData. Late-stage feature (modalities processed independently first). Patient-level metadata propagation across modalities. |
| **Cross-modal queries** | `sct query --project ibd_spatial --patient P001 --modalities rna,visium` returns available data per patient across modalities. Enables multi-omic analysis planning. | Medium | Requires patient metadata attachment (Phase 2 output) and assembly. Returns JSON inventory. |
| **Command chaining with pipe support** | `sct qc run --input raw.h5ad | sct preprocess run` -- commands accept input from previous command's JSON output (specifically the `output_file` field). | Medium | Each command outputs `{"output_file": "...", "status": "success", ...}`. Next command reads `--input` from stdin JSON. |
| **GPU/resource negotiation** | `sct preprocess run --gpu auto` detects GPU availability, selects rapids-singlecell vs scanpy, reports which backend was used. Existing `_gpu.py` does detection but doesn't expose it to agents. | Low | Wraps existing GPU detection. Return `{"backend": "rapids", "gpu": "A100", "memory_gb": 40}` in output metadata. |
| **Skill file for agent guidance** | Ship `.sct-skills.md` with structured YAML frontmatter documenting agent-specific best practices: "Always run `sct validate` after `sct preprocess`", "Use `--subsample` for datasets >1M cells". | Low | Static file, but high value. Agents read this once and follow conventions. Like CLI-Anything's skill files. |

## Anti-Features

Features to explicitly NOT build. Each represents a tempting trap.

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| **Interactive TUI/wizard mode** | Agents cannot use interactive interfaces. Humans who want interactivity use Jupyter notebooks (the existing workflow). A TUI serves neither audience well. | Comprehensive `--help` and `sct describe` for discoverability. |
| **Full workflow engine / Snakemake replacement** | PROJECT.md explicitly scopes this out. sc_tools CLI complements Snakemake, not replaces it. Building a workflow engine is a multi-year project that competes with Nextflow, Snakemake, WDL. | Composable single commands. Agents orchestrate the workflow; the CLI provides the building blocks. |
| **Web dashboard / GUI** | Out of scope (PROJECT.md). HTML reports already exist. A dashboard adds auth, state management, frontend framework -- massive scope. | `sct report generate` produces static HTML reports. Agents serve these to humans. |
| **Database-first provenance** | Schema is unstable (PROJECT.md: "frequent migrations, weak typing"). Committing to a DB schema now means constant migrations. | JSON sidecar files. Migrate to DB once model stabilizes from actual usage patterns. |
| **Custom data format** | AnnData/MuData are the scverse standard. Inventing a custom format fragments the ecosystem and forces conversion. CellRanger's proprietary formats are a cautionary tale. | AnnData in, AnnData out. MuData for multi-modal. Standard formats agents and humans both know. |
| **Automatic parameter optimization** | Tempting to auto-tune scVI latent dims, Leiden resolution, etc. But this hides decisions from the agent/researcher, makes results non-reproducible, and is a research problem not an engineering one. | Expose parameters explicitly. Agents + researchers make decisions. Provide sensible defaults documented in `--help`. |
| **Real-time streaming progress bars** | Rich progress bars are human conveniences that pollute agent stdout. ANSI escape codes in JSON output break parsing. | `--progress` flag writes progress JSON to stderr (optional). Silent by default for agents. |
| **Plugin/extension system** | Premature architecture. The integration methods, cell typing methods, and platforms are already pluggable in sc_tools internals. A formal plugin API adds maintenance burden without current need. | Use sc_tools' existing internal extensibility. Add new methods to the library, CLI wraps them automatically. |

## Feature Dependencies

```
Structured JSON output ──────────────────────────┐
                                                   │
Self-describing help ──────────────────────────────┤
                                                   ├──► ALL other features depend on these
Semantic exit codes ───────────────────────────────┤
                                                   │
Structured error reporting ────────────────────────┘

Checkpoint validation ──► Phase-aware pipeline status
                                │
                                ▼
                    JSON sidecar provenance (enriches status with lineage)

Memory estimation ──► Dry-run mode (estimation is a subset of dry-run)

Pre-computed benchmark ingestion ──► (independent, but benefits from provenance)

Patient metadata attachment ──► MuData assembly ──► Cross-modal queries

GPU/resource negotiation ──► Memory estimation (GPU changes memory profile)
```

## MVP Recommendation

### Phase 1: Agent Interface Foundation
Build the features that make the CLI usable by agents at all:

1. **Structured JSON output** -- without this, nothing else matters
2. **Self-describing help** -- agents need discovery
3. **Semantic exit codes + structured error reporting** -- agents need to handle failures
4. **Non-interactive execution** -- agents can't use prompts
5. **Verb-noun command structure** -- stable command grammar
6. **Phase-aware pipeline status** -- wraps existing registry, immediate value

### Phase 2: Reliability and Safety
Features that make agent workflows robust:

7. **Checkpoint validation** -- verify before proceeding
8. **Idempotent operations** -- safe retries
9. **Dry-run mode** -- planning before execution
10. **Memory estimation** -- prevent OOM on large datasets
11. **GPU/resource negotiation** -- expose existing capability

### Phase 3: Provenance and Benchmarking
The domain-specific differentiators:

12. **JSON sidecar provenance** -- reproducibility
13. **Pre-computed benchmark ingestion** -- solves the 44G OOM problem (core motivation)
14. **Skill file** -- agent guidance

### Phase 4: Multi-Modal
Late-stage, requires modalities to be independently processed first:

15. **Patient metadata attachment** (wraps existing Phase 2 capabilities)
16. **MuData assembly** -- link modalities
17. **Cross-modal queries** -- multi-omic analysis planning

**Defer indefinitely:** Command chaining with pipes. Nice to have, but agents can read JSON output and construct the next command themselves. The complexity of stdin/stdout coordination isn't worth it when the agent is the orchestrator.

## Sources

- [CLI-Anything: Making ALL Software Agent-Native](https://github.com/HKUDS/CLI-Anything) -- agent-native CLI generation patterns
- [You Need to Rewrite Your CLI for AI Agents](https://justin.poehnelt.com/posts/rewrite-your-cli-for-ai-agents/) -- structured output, schema introspection, skill files
- [Keep the Terminal Relevant: Patterns for AI Agent Driven CLIs (InfoQ)](https://www.infoq.com/articles/ai-agent-cli/) -- exit codes, dry-run, MCP exposure, output contracts
- [Nextflow Data Lineage Tutorial](https://nextflow.io/docs/latest/tutorials/data-lineage.html) -- provenance tracking patterns in scientific pipelines
- [Nextflow Feature Highlights: Data Lineage (Seqera)](https://seqera.io/blog/nextflow-updates-strict-syntax-data-lineage/) -- Nextflow 25.04 native lineage
- [scanpy-scripts (PyPI)](https://pypi.org/project/scanpy-scripts/) -- existing scanpy CLI wrapper (limited, no structured output)
- [MuData documentation](https://mudata.readthedocs.io/en/latest/io/mudata.html) -- multi-modal data container architecture
- [muon documentation](https://muon.readthedocs.io/en/latest/io/mudata.html) -- MuData annotation management
- [gh CLI Manual](https://cli.github.com/manual/gh_help_reference) -- reference implementation for `--json` flag pattern
