# Requirements — v2.0 Report Plots & Sample Concat

## v2.0 Requirements

### Maintenance (Prerequisites)

- [ ] **MAINT-01**: `plotly` promoted from optional `[pipeline]` dependency to base `dependencies` in `pyproject.toml`
- [ ] **MAINT-02**: Stale CDN pin (`plotly-2.27.0.min.js`) replaced with `plotly-3.4.0.min.js` (pinned version matching installed plotly.py 6.6.0) in `report_utils.py`

### Sample Concatenation

- [ ] **CONCAT-01**: User can run `sct concat` with a list of same-modality h5ad paths to produce a merged h5ad
- [ ] **CONCAT-02**: `sct concat` preserves `uns["spatial"]` across all input samples via `uns_merge="unique"` so downstream spatial plots render correctly
- [ ] **CONCAT-03**: `sct concat` writes a `.provenance.json` sidecar recording input files with SHA256 checksums, sc_tools version, and timestamp
- [ ] **CONCAT-04**: `sct concat` is registered as a pipeline phase (`concat`) in `pipeline.py` between `ingest_load` and `qc_filter`, visible in `sct status` output

### QC Report Spatial Plots

- [ ] **QC-01**: `sct qc run` computes and stores `tpm_worthy` (boolean) in `adata.obs` — a spot/cell is TPM-worthy if it meets count and gene thresholds sufficient for TPM normalization
- [ ] **QC-02**: `sct qc run` computes a per-sample QC score: passing samples (0–100) = 0.4 × (% TPM-worthy spots) + 0.3 × min(median_genes_in_TPM_spots / 5000, 1) × 100 + 0.3 × MT penalty (0% MT = full credit, ≥20% MT = 0); failed samples receive a negative score by failure mode: complete failure = −100, RNA degradation = −60, insufficient area = −40, sparse tissue = −30, low transcription = −20
- [ ] **QC-03**: `sct qc run` flags and excludes samples with fewer than 100 TPM-worthy spots as unusable, reporting them separately in the CLIResult
- [ ] **QC-04**: The pre-filter QC report renders 4 spatial panels per sample: `log1p_counts`, `log1p_genes`, `pct_counts_mt`, and `tpm_worthy` (TPM-worthy spot overlay), with samples ordered by ascending QC score (worst first)
- [ ] **QC-05**: Spatial panels are skipped gracefully for non-spatial modalities (no `obsm["spatial"]`) with a warning rather than an error

### Integration Report UMAPs

- [ ] **INT-01**: The post-integration report includes `subject_id`/patient as a UMAP color key alongside existing `leiden`, `batch`, and `sample` keys
- [ ] **INT-02**: The integration report UMAP grid uses a bio-vs-batch split layout: batch-correction panels (batch, sample) on one side, biology-conservation panels (leiden, celltype if present, subject_id) on the other
- [ ] **INT-03**: Bio obs columns used for scib metrics are auto-detected from the provenance sidecar (`.provenance.json`) rather than requiring explicit CLI flags

### Celltype Report Spatial Plots

- [ ] **CELL-01**: The post-celltyping report renders one aggregate spatial plot per sample colored by cell type assignment
- [ ] **CELL-02**: The post-celltyping report renders individual per-cell-type spatial overlays for each cell type, with a configurable `--max-celltypes` cap (default: 20) to prevent HTML size blowup
- [ ] **CELL-03**: Spatial plots in the celltype report are skipped gracefully for non-spatial modalities with a warning

## Future Requirements

- `sct estimate concat` — memory pre-check before running concat on large samples (deferred, assess after v2.0 usage)
- Backed-mode concat for VisiumHD samples at scale (deferred to v2.1)
- Per-spot QC threshold display in report footer (deferred)
- Plotly CDN offline fallback for HPC environments without internet (deferred)

## Out of Scope

- Cross-modality concat (e.g., Visium + IMC in one call) — same-modality only in v2.0
- Batch pipeline / manifest-driven `sct run` across multiple samples — identified as premature
- Interactive spatial plots (Plotly/Bokeh) — base64 PNG embedding is the established pattern; interactive spatial at 2.5M cells produces unacceptably large HTML
- New integration methods — use existing sc_tools implementations

## Traceability

| REQ-ID | Phase | Status |
|--------|-------|--------|
| MAINT-01 | Phase 9 | Pending |
| MAINT-02 | Phase 9 | Pending |
| CONCAT-01 | Phase 9 | Pending |
| CONCAT-02 | Phase 9 | Pending |
| CONCAT-03 | Phase 9 | Pending |
| CONCAT-04 | Phase 9 | Pending |
| QC-01 | Phase 10 | Pending |
| QC-02 | Phase 10 | Pending |
| QC-03 | Phase 10 | Pending |
| QC-04 | Phase 10 | Pending |
| QC-05 | Phase 10 | Pending |
| CELL-01 | Phase 11 | Pending |
| CELL-02 | Phase 11 | Pending |
| CELL-03 | Phase 11 | Pending |
| INT-01 | Phase 12 | Pending |
| INT-02 | Phase 12 | Pending |
| INT-03 | Phase 12 | Pending |
