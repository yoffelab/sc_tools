# Roadmap: sc_tools Agent-Native CLI

## Milestones

- ✅ **v1.0 Agent-Native CLI** - Phases 1-8 (shipped 2026-03-24)
- 🚧 **v2.0 Report Plots & Sample Concat** - Phases 9-12 (in progress)

## Phases

<details>
<summary>v1.0 Agent-Native CLI (Phases 1-8) - SHIPPED 2026-03-24</summary>

- [x] **Phase 1: sc_tools Internals Fixes** - Fix integration benchmark bugs, subsampling bias, and recipe ordering
- [x] **Phase 2: CLI Foundation** - Typer-based sct entry point with CLIResult envelope and structured error reporting
- [x] **Phase 3: CLI Core Commands** - sct qc run, preprocess run, benchmark integration, validate, status
- [x] **Phase 4: CLI Self-Discovery** - list-commands, describe, schema for machine-readable introspection
- [x] **Phase 5: Provenance & Reproducibility** - JSON sidecar provenance files with lineage tracing
- [x] **Phase 6: Scientific Gaps** - Pseudobulk DE, marker validation, subject metadata, panel-aware cell typing
- [x] **Phase 7: Memory Safety** - IO Gateway with tiered loading, sct estimate, dry-run mode
- [x] **Phase 8: Multi-Omic Assembly** - MuData assembly, patient metadata join, cross-modal queries

</details>

### 🚧 v2.0 Report Plots & Sample Concat (In Progress)

**Milestone Goal:** Add spatial and UMAP plots to QC/integration/celltype HTML reports, and expose sample concatenation as a first-class CLI command.

- [ ] **Phase 9: Sample Concatenation & Maintenance** - sct concat CLI command with pipeline registration, plus plotly dependency fixes
- [ ] **Phase 10: QC Spatial Plots** - Per-sample spatial panels in pre-filter QC report with TPM-worthy scoring and sample ranking
- [ ] **Phase 11: Celltype Spatial Plots** - Per-sample and per-celltype spatial overlays in post-celltyping report
- [ ] **Phase 12: Integration Report UMAPs** - Subject-level UMAP coloring with bio-vs-batch split layout and auto-detected obs columns

## Phase Details

### Phase 9: Sample Concatenation & Maintenance
**Goal**: Users can merge same-modality samples via CLI and all report infrastructure dependencies are resolved
**Depends on**: Phase 8 (v1.0 complete)
**Requirements**: MAINT-01, MAINT-02, CONCAT-01, CONCAT-02, CONCAT-03, CONCAT-04
**Success Criteria** (what must be TRUE):
  1. User can run `sct concat --input s1.h5ad s2.h5ad --output merged.h5ad` and get a valid merged AnnData with all spatial coordinates preserved
  2. `sct concat` output includes a `.provenance.json` sidecar with SHA256 checksums of all input files
  3. `sct status` shows `concat` as a recognized pipeline phase positioned between `ingest_load` and `qc_filter`
  4. `import plotly` succeeds in a base sc_tools install (no `[pipeline]` extras needed) and report HTML loads current Plotly JS
**Plans**: TBD

### Phase 10: QC Spatial Plots
**Goal**: Pre-filter QC report shows spatial distribution of key metrics per sample, ordered by QC quality, enabling users to visually identify problematic samples before filtering
**Depends on**: Phase 9
**Requirements**: QC-01, QC-02, QC-03, QC-04, QC-05
**Success Criteria** (what must be TRUE):
  1. After `sct qc run`, the output AnnData contains `tpm_worthy` in obs and users can inspect which spots meet TPM normalization thresholds
  2. `sct qc run` CLIResult reports a per-sample QC score (0-100 for passing, negative for failed) and separately lists samples excluded for having fewer than 100 TPM-worthy spots
  3. The pre-filter QC HTML report displays 4 spatial panels per sample (log1p_counts, log1p_genes, pct_counts_mt, tpm_worthy) with worst-scoring samples shown first
  4. For non-spatial modalities (no `obsm["spatial"]`), the report generates successfully with spatial sections skipped and a warning logged
**Plans**: TBD
**UI hint**: yes

### Phase 11: Celltype Spatial Plots
**Goal**: Post-celltyping report shows spatial distribution of cell type assignments per sample, enabling users to assess whether cell type calls are spatially coherent
**Depends on**: Phase 10
**Requirements**: CELL-01, CELL-02, CELL-03
**Success Criteria** (what must be TRUE):
  1. The post-celltyping HTML report renders one aggregate spatial plot per sample colored by cell type assignment with a consistent color palette across all samples
  2. The post-celltyping report renders individual per-cell-type spatial overlays, capped at `--max-celltypes` (default 20) to prevent HTML bloat
  3. For non-spatial modalities, the celltype report generates successfully with spatial sections skipped and a warning logged
**Plans**: TBD
**UI hint**: yes

### Phase 12: Integration Report UMAPs
**Goal**: Post-integration report shows UMAP embeddings colored by both batch-correction and biology-conservation dimensions, enabling users to visually assess integration quality
**Depends on**: Phase 9
**Requirements**: INT-01, INT-02, INT-03
**Success Criteria** (what must be TRUE):
  1. The post-integration HTML report includes `subject_id`/patient as a UMAP color key when the column exists in the AnnData
  2. The UMAP grid uses a bio-vs-batch split layout: batch-correction panels (batch, sample) separated from biology-conservation panels (leiden, celltype, subject_id)
  3. Bio obs columns used for coloring are auto-detected from the `.provenance.json` sidecar without requiring explicit CLI flags
**Plans**: TBD
**UI hint**: yes

## Progress

**Execution Order:**
Phases execute in numeric order: 9 → 10 → 11 → 12
Note: Phase 12 depends on Phase 9 (not 11), so it could execute in parallel with 10-11 if needed.

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 1. Internals Fixes | v1.0 | 2/2 | Complete | 2026-03-21 |
| 2. CLI Foundation | v1.0 | 2/2 | Complete | 2026-03-21 |
| 3. CLI Core Commands | v1.0 | 3/3 | Complete | 2026-03-22 |
| 4. CLI Self-Discovery | v1.0 | 1/1 | Complete | 2026-03-23 |
| 5. Provenance | v1.0 | 2/2 | Complete | 2026-03-23 |
| 6. Scientific Gaps | v1.0 | 4/4 | Complete | 2026-03-24 |
| 7. Memory Safety | v1.0 | 2/2 | Complete | 2026-03-24 |
| 8. Multi-Omic Assembly | v1.0 | 2/2 | Complete | 2026-03-24 |
| 9. Sample Concat & Maintenance | v2.0 | 0/? | Not started | - |
| 10. QC Spatial Plots | v2.0 | 0/? | Not started | - |
| 11. Celltype Spatial Plots | v2.0 | 0/? | Not started | - |
| 12. Integration Report UMAPs | v2.0 | 0/? | Not started | - |
