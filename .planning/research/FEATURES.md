# Feature Research

**Domain:** Spatial/UMAP plot-rich HTML reports + sample concatenation CLI for sc_tools v2.0
**Researched:** 2026-03-27
**Confidence:** HIGH (based on direct codebase inspection + domain knowledge of scanpy/squidpy conventions)

---

## Current State Baseline (What Already Exists)

Before listing features to build, the existing report infrastructure must be understood to avoid re-implementing what is already there.

### Pre-filter QC report (`generate_pre_filter_report`)
Already has: QC 2x2 histogram grid, counts-vs-genes scatter, violin distributions, per-sample violin grouped (linear + log10), cross-sample comparison bars (linear + log10), QC scatter matrix, `%MT` per sample. Template has a `spatial_multipage` slot but **the generation function never populates it** — it is dead code. No spatial plots are generated for pre-filter.

### Post-filter QC report (`generate_post_filter_report`)
Already has: all pre-filter plots plus HVG/SVG feature selection plots, AND a `spatial_counts_plots` list of per-sample `log1p(total_counts)` spatial overlays (one per sample, via `sc.pl.spatial`). This is the **only existing spatial plot** in the report pipeline, and only for post-filter. No `log1p_genes`, no `%mt`, no pass/fail spatial overlay.

### Post-integration QC report (`generate_post_integration_report`)
Already has: a `umap_grid` (single composite PNG colored by sample/batch/cluster/celltype), per-embedding UMAP comparison grid (one panel per integration method, colored by sample only), cluster distribution stacked bar, integration metrics radar (Plotly), batch-vs-bio scatter (Plotly). Missing: UMAP panels for `subject_id`/patient are not auto-detected. The single `umap_grid` bundles all colorings into one figure without batch-vs-bio structural separation.

### Post-celltyping QC report (`generate_post_celltyping_report`)
Already has: celltype abundance bar per sample, UMAP grid (celltype/celltype_broad/leiden/batch/sample), cluster distribution by celltype, marker dotplot, marker validation table. **Missing: all spatial plots** — no per-cell-type or per-sample spatial overlays of assigned cell types.

### Concat
No `sct concat` CLI command exists. No backend concat wrapper in `sc_tools` beyond what `anndata.concat` provides natively. The pipeline DAG in `pipeline.py` has no concat phase registered.

---

## Feature Landscape by Report Category

---

### A. QC Report — Spatial Plots (Pre-filter)

#### Table Stakes (Users Expect These)

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| Spatial plot: `log1p_total_counts` per sample, pre-filter | First thing a Visium analyst checks — are counts distributed uniformly or is there tissue-edge dropout? | LOW | `sc.pl.spatial` already called in post-filter; replicate for pre-filter in `generate_pre_filter_report`. Uses the dead `spatial_multipage` template slot. |
| Spatial plot: `log1p_n_genes_by_counts` per sample, pre-filter | Gene detection spatial pattern reveals tissue quality distinct from count depth | LOW | Requires computing `log1p_n_genes_by_counts` on obs if not present, then `sc.pl.spatial`. Column already computable from `n_genes_by_counts`. |
| Spatial plot: `pct_counts_mt` per sample, pre-filter | MT contamination is often spatially structured (necrotic core, tissue edges); must see spatial pattern not just distribution | LOW | Column already computed by `calculate_qc_metrics`; only the spatial visualization is missing. |
| Spatial plot: pass/fail per spot, pre-filter | Critical for QC review — shows which spots were removed and whether the spatial pattern of removal is sensible (edge effect vs. random noise) | MEDIUM | Requires adding a binary `qc_pass_spot` obs column during `sct qc run`. Per-spot pass/fail is not currently stored; `classify_samples` operates at sample level only. |
| Sample ordering by QC score in spatial section | Failed samples should appear first (worst first) so reviewers immediately see the problem cases | LOW | Currently the sample loop uses `sorted()` alphabetically. Change to sort by `classified` QC pass/score before the loop. |

#### Differentiators

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| 4-panel spatial grid per sample (counts, genes, mt, pass/fail) | Single-glance spatial quality per sample — matches how analysts actually review Visium data | MEDIUM | Generate a 2x2 subplot per sample using `plt.subplots(2,2)` with `sc.pl.spatial(..., ax=ax)` in each panel. All 4 metrics in one figure per sample rather than 4 separate loops. |
| Common color scale across samples for each metric | Enables direct visual comparison across samples — without common scale, a "good-looking" sample may just have a compressed colormap | MEDIUM | Compute cohort-level 2nd–98th percentile vmin/vmax for each continuous metric before the sample loop. Pass `vmin`/`vmax` to each `sc.pl.spatial` call. |

#### Anti-Features

| Feature | Why Requested | Why Problematic | Alternative |
|---------|---------------|-----------------|-------------|
| Spatial plots for non-Visium modalities (IMC, Xenium) in this report | IMC/Xenium also have spatial coordinates | `sc.pl.spatial` requires `adata.uns['spatial'][library_id]['images']`; this only exists for Visium/VisiumHD. IMC uses `X_spatial` in obsm without tissue images. Calling `sc.pl.spatial` on IMC data raises or silently produces garbage. | Guard all QC spatial sections with `"spatial" in adata.uns and library_id in adata.uns.get("spatial", {})`. Skip silently for non-Visium modalities. IMC/Xenium spatial QC is a separate future feature. |
| Recomputing per-spot pass/fail thresholds inside report generation | "Just show which spots fail the QC thresholds" | Re-implementing the filter logic in the report creates a divergence risk — the report may show different pass/fail than what was actually applied in `sct qc run`. | Store `qc_pass_spot` in `adata.obs` during `sct qc run` before the filter step. Report reads it from obs. No recomputation in report generation. |
| Embedding individual spatial plot PNGs as separate files | "Easier to share individual plots" | Self-contained HTML is the design contract for sc_tools reports. Per-file outputs fragment the report. | Base64 inline. For very large cohorts (>20 samples), produce a grid of thumbnails in one figure rather than full-resolution per-sample figures. |

---

### B. Integration Report — UMAP Plots

#### Table Stakes (Users Expect These)

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| UMAP colored by Leiden clusters | Leiden is the primary cluster label used downstream; cluster structure must be visible in the integration report | LOW | `leiden` key already in `color_keys` list if present in `adata.obs`. Already in `umap_grid`. Gap: it is bundled with all other colorings in one multi-panel figure with no structural separation. |
| UMAP colored by batch key | Batch mixing is the primary integration QA question; a well-integrated dataset should show batch labels mixed across the UMAP | LOW | `batch_key` already in `color_keys`. Structural separation from bio colorings is the missing piece. |
| UMAP colored by `library_id` / sample | Sample-level mixing — distinct from batch when multiple libraries share a batch label | LOW | `sample_col` already in `color_keys`. |
| UMAP colored by `subject_id` / patient | Subject-level biological signal — critical for multi-patient cohorts to distinguish patient effect from batch effect | MEDIUM | `subject_id` is the v1.0 subject-level metadata model field (Phase 6). Must be added to `color_keys` auto-detection in `generate_post_integration_report`. Currently not auto-detected. |
| UMAP colored by bio obs columns used in benchmark metrics | The bio metrics (ARI/NMI/ASW-cell) were computed against `celltype_key` and `celltype_broad`; UMAP must show those labels for the metrics to be interpretable | LOW | `celltype_key` already conditionally added when `has_celltype=True`. `celltype_broad` already added if present. This is table stakes only if celltype annotation has been run before integration report generation (which is the post-celltyping report scenario). |

#### Differentiators

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Split UMAP layout: batch-side vs bio-side | Two columns — left: batch/sample/subject UMAPs (batch QC), right: leiden/celltype UMAPs (bio QC) — makes the report self-interpreting without reading captions | MEDIUM | Currently all color keys concatenated into one `qc_umap_grid` figure. Requires splitting `color_keys` into two lists and either two figures or a structured subplot layout with a dividing header. |
| Per-embedding UMAP colored by Leiden (not just sample) | Current per-embedding grid is colored only by `sample`; coloring by Leiden shows whether clusters are consistent across integration methods, not just whether samples mix | MEDIUM | Requires a second `qc_embedding_umap_grid` call with `color_key=cluster_key`. Doubles compute time for the grid. Run after confirming cluster key is available. |

#### Anti-Features

| Feature | Why Requested | Why Problematic | Alternative |
|---------|---------------|-----------------|-------------|
| Computing UMAP fresh during report generation | "Show the best UMAP" | Report generation should be read-only and fast. Recomputing UMAP at 2.5M cells (VisiumHD) takes 30+ minutes and changes analysis state. | Require `X_umap` to be present in `adata.obsm` before generating integration report. Emit a warning and skip UMAP section if absent. UMAP computation belongs in `sct preprocess run`. |
| One UMAP panel per integration method per color key | "Comprehensive — show everything" | 5 methods × 4 color keys = 20 UMAP panels; unreadable in HTML; each UMAP computation for a new embedding can take minutes | Per-embedding grid shows one color key (sample); main UMAP grid shows all color keys for the selected/best method only. |

---

### C. Celltype Report — Spatial Plots

#### Table Stakes (Users Expect These)

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| Spatial overlay of `celltype` label per sample | The fundamental output of cell type annotation for spatial data — seeing where cell types are spatially distributed is the primary value of spatial transcriptomics | HIGH | Requires `plot_spatial_categorical` (already in `sc_tools/pl/spatial.py`) called per sample, looping over `library_id` values. Result: list of per-sample figures, each showing categorical cell type overlay on tissue. Legend crowding with many cell types requires a pre-built palette and `frameon=False`. |
| Consistent color palette across samples | Cell types must have the same color across all sample panels to be comparable | LOW | Pre-compute a `palette` dict from all unique `celltype` values in `adata.obs[celltype_key]` before the sample loop. Pass `palette=` to each `plot_spatial_categorical` call. |
| Modality guard — Visium/VisiumHD only | Spatial cell type overlays only work where tissue image + spot coordinates exist in `uns['spatial']` | LOW | Check `"spatial" in adata.uns` and modality is in `{"visium", "visium_hd", "visium_hd_cell"}` before generating spatial section. For non-spatial modalities emit an info log and skip section. |

#### Differentiators

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| `celltype_broad` spatial overlay if present | Broad categories (immune/stromal/tumor) easier to interpret at a glance than fine-grained types | LOW | Same pattern as `celltype`; add a second spatial section with `celltype_broad` key if the column exists in `adata.obs`. |
| Spatial overlay of `leiden` cluster per sample | Shows concordance between Leiden clusters and annotated cell types in tissue space; mismatches reveal annotation uncertainty | MEDIUM | Same loop pattern as celltype; `plot_spatial_categorical` with `color=cluster_key`. Useful for QA of annotation quality. |

#### Anti-Features

| Feature | Why Requested | Why Problematic | Alternative |
|---------|---------------|-----------------|-------------|
| Per-cell-type spatial overlays, all cell types inline | "Show where each cell type is distributed" — high scientific value | N_samples=10 × N_celltypes=15 = 150 plots; at 150 DPI each PNG is ~1MB; HTML grows to 150+ MB; browsers hang or crash. | MVP: aggregate categorical overlay only (all cell types in one plot per sample). Post-MVP: accordion/collapsible per-cell-type section limited to top 8 cell types by abundance. Per-cell-type PDF as a separate artifact is a better format for this content. |
| Spatial plots for IMC without image backend | IMC has spatial coordinates | `sc.pl.spatial` requires `uns['spatial'][library_id]['images']`; IMC uses channel images, not H&E. `sc.pl.spatial` will fail. `pl.spatial.plot_imc_composite` (already in the codebase) is the correct call for IMC but takes different arguments. | Guard report spatial section with modality check. IMC cell type spatial visualization uses a separate code path and is a future feature. |

---

### D. `sct concat` Command

#### Table Stakes (Users Expect These)

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| Concatenate multiple same-modality h5ad files into one | Multi-sample QC and integration require a single concatenated AnnData; this is the entry point for every multi-sample pipeline | MEDIUM | Wrap `anndata.concat(adatas, label=batch_key_col, keys=sample_ids, uns_merge="unique")`. Accept a list of paths via `--input` (repeatable) or `--from-dir` with a glob pattern. |
| `obs_names_make_unique` / obs index deduplication | After concat, obs indices from different samples collide; scanpy and scVI require unique obs names | LOW | Use `index_unique="-"` in `anndata.concat` to suffix obs names with the batch label. Must be on by default; document that this changes obs index format. |
| Batch key injection into `adata.obs` | Integration methods (scVI, Harmony, BBKNN) require a `batch_key` obs column; concat must create it | LOW | `anndata.concat(label=batch_key_col, keys=sample_ids)` automatically creates the label column. Configurable via `--batch-key library_id` (default). |
| `uns['spatial']` preservation across all input samples | `uns['spatial']` contains per-sample H&E image and scalefactor data; losing it silently breaks all downstream spatial plots | MEDIUM | Use `uns_merge="unique"` in `anndata.concat`. Verify post-concat that all input sample IDs appear as keys in `output.uns['spatial']`. Fail with a clear error if any sample's spatial metadata is lost. This is the correctness-critical requirement for Visium data. |
| Layer preservation (`layers['counts']`) | Raw counts in `layers['counts']` are the standard sc_tools convention; scVI requires raw counts and reads them from this layer | LOW | `anndata.concat` preserves layers by default when all inputs have the same layer set. Emit a warning if layer sets are mismatched across samples (e.g., one sample has `layers['counts']`, another does not). |
| CLIResult JSON stdout with structured data | Agent-native output consistent with all other `sct` commands | LOW | Standard CLIResult pattern. `data` dict must include: `n_obs_total`, `n_vars`, `n_samples`, `samples` (list of library IDs), `output_path`. |
| Provenance sidecar (`.provenance.json`) | Phase 5 requirement — all `sct` commands write provenance | LOW | Record: input file paths + checksums, `--batch-key` value, `uns_merge` strategy, `index_unique` setting, sc_tools version, timestamp, runtime, peak memory. Same pattern as `sct qc run`. |
| Pipeline phase registration in `pipeline.py` | `sct status` must show concat phase; downstream commands (`sct qc run`) must list it as a prerequisite for multi-sample data | MEDIUM | Add `PhaseSpec` for `concat` in `pipeline.py` STANDARD_PHASES. Checkpoint: output h5ad exists. Required obs after concat: `library_id` (or configured `batch_key`). Required uns for Visium: `spatial` with all sample keys. |

#### Differentiators

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Pre-flight validation: per-sample obs counts and var space consistency | Catch data problems before concat (empty samples, mismatched gene sets) with clear error messages rather than silent corruption | LOW | Before `anndata.concat`: check `adata.n_obs > 0` for each input; check all inputs have identical `adata.var_names` (or warn on mismatch); log per-sample obs count in the CLIResult `data.samples` list. |
| Modality consistency validation | Catch modality mismatch early — mixing Visium and VisiumHD var names produces an invalid AnnData that breaks all downstream tools | LOW | Check that all input h5ad have the same `adata.uns.get('sc_tools_modality')` value. Warn on mismatch; add `--strict-modality` flag to fail hard. |
| Lexicographic sort guarantee for `--from-dir` input discovery | Reproducible sample ordering when discovering files from a directory; order affects batch key index assignment | LOW | Sort glob results lexicographically before concat. Document that alphabetical order is the stable default. |

#### Anti-Features

| Feature | Why Requested | Why Problematic | Alternative |
|---------|---------------|-----------------|-------------|
| Cross-modality concat (Visium + IMC in one `sct concat`) | "One command for all modalities" | Different var spaces, different spatial coordinate conventions, different obs structures — concat produces a semantically invalid AnnData that breaks every downstream tool | Enforce same-modality constraint (same `sc_tools_modality` tag in uns). Cross-modality assembly is MuData via `sct assemble` (Phase 8, already scoped). |
| Selecting integration method or running integration during concat | "Make it smart — concat and integrate in one step" | Concat is a pure data operation. Coupling it to integration creates hidden state and makes the command non-idempotent. | Keep concat a data-only command. Integration method selection happens in `sct preprocess run` or `sct benchmark integration`. |
| Storing concat output in registry as a phase output | "Provenance tracking" | Registry schema has 7 unapplied migrations as of v1.0; adding a new phase to the registry before the schema stabilizes will create another migration | File-based provenance sidecar only. Registry registration can be added as a migration after the schema stabilizes. |
| Memory-efficient backed-mode concat for VisiumHD at MVP | "VisiumHD files are 11–25G each" | Backed-mode concat with `anndata.concat` requires all inputs to be in backed mode, which has different API constraints and limited write support. Implementing this correctly is HIGH complexity. | For MVP, document the memory requirement (N_samples × ~4× file size) and let the IO Gateway (`sct estimate`) warn before exec. Backed-mode concat is a v1.x feature. |

---

## Feature Dependencies

```
sct concat
    └──enables──> sct qc run (multi-sample, pre-filter)
                      └──enables──> generate_pre_filter_report
                                        └──needs──> spatial plots (currently dead slot)
                                        └──needs──> pass/fail spatial overlay
                                                        └──requires──> qc_pass_spot in adata.obs
                                                                           └──requires──> sct qc run adds per-spot column before filter

generate_post_integration_report
    └──needs──> UMAP colored by subject_id
                    └──requires──> subject_id in adata.obs
                                       └──satisfied by v1.0 Phase 6 subject metadata model
    └──needs──> split batch-vs-bio UMAP layout
                    └──requires──> refactored qc_umap_grid or two separate figure calls

generate_post_celltyping_report
    └──needs──> spatial per-sample celltype overlay
                    └──requires──> plot_spatial_categorical (already in sc_tools/pl/spatial.py)
                    └──requires──> spatial image data in uns['spatial'] (modality guard needed)
                    └──conflicts with HTML size if N_samples x N_celltypes > 50 plots
```

### Dependency Notes

- **`qc_pass_spot` must be produced by `sct qc run`, not re-derived in the report.** The sample-level `classify_samples` DataFrame does not have per-spot decisions. The per-spot mask must be stored in `adata.obs['qc_pass_spot']` during the QC phase. Report reads it from obs; no filtering logic in report code.

- **Spatial plots require `uns['spatial']` preserved through concat.** If `sct concat` does not use `uns_merge="unique"`, all spatial plot generation silently fails because `adata.uns['spatial']` will be missing or partial after concat. Concat correctness is a prerequisite for all QC spatial features.

- **Integration UMAP panels depend on `subject_id` being present in obs.** The `subject_id` field was added in Phase 6 (validated 2026-03-24). It is only present in samples processed with the v1.0 subject metadata model. The integration report should auto-detect it (`if "subject_id" in adata.obs.columns`) rather than require it.

- **Celltype spatial section must guard on modality.** `sc.pl.spatial` requires `adata.uns['spatial'][library_id]['images']`. This only exists for Visium and VisiumHD. All spatial plot generation must be gated on modality via `get_modality_terms(modality)` which already tracks `_CELL_MODALITIES`.

---

## MVP Definition for v2.0

### Launch With (v1 of this milestone)

QC spatial plots (pre-filter report):
- [ ] `log1p_total_counts`, `log1p_n_genes_by_counts`, `pct_counts_mt` spatial overlays per sample — replaces the dead `spatial_multipage` slot
- [ ] Sample ordering by QC score (failed samples first) in spatial section
- [ ] `qc_pass_spot` obs column added during `sct qc run` before filter step

Integration UMAP:
- [ ] Add `subject_id` to `color_keys` auto-detection in `generate_post_integration_report`
- [ ] Pass/fail spatial overlay reads `qc_pass_spot` from obs when present

Celltype spatial:
- [ ] Aggregate per-sample spatial overlay of `celltype` categorical label — one figure per sample, consistent palette
- [ ] Modality guard (Visium/VisiumHD only)

`sct concat`:
- [ ] `sct concat --input <file> [<file>...] --output <path> --batch-key library_id`
- [ ] obs dedup (`index_unique="-"`), uns spatial merge (`uns_merge="unique"`), layer passthrough
- [ ] Pre-flight: per-sample obs count validation, var space consistency check
- [ ] CLIResult JSON with `n_obs_total`, `n_vars`, `n_samples`, `samples`, `output_path`
- [ ] Provenance sidecar
- [ ] `PhaseSpec` registration in `pipeline.py`

### Add After Validation (v1.x)

- [ ] Common color scale across samples for spatial metrics (cohort-level vmin/vmax)
- [ ] 4-panel spatial grid per sample (2x2: counts/genes/mt/pass-fail) in pre-filter report
- [ ] Split batch-vs-bio UMAP layout in integration report (two-column figure)
- [ ] Per-celltype spatial overlays in celltype report — accordion, top 8 by abundance only
- [ ] `celltype_broad` spatial overlay in celltype report
- [ ] `sct concat --from-dir` with modality validation and lexicographic sort

### Future Consideration (v2+)

- [ ] Per-embedding UMAP colored by Leiden (doubles grid compute time; validate demand first)
- [ ] Spatial plots for IMC (requires `plot_imc_composite` integration into reports; different code path)
- [ ] Spatial plots for Xenium/CosMx (different coordinate conventions; needs separate spatial rendering path)
- [ ] Memory-efficient backed-mode concat for VisiumHD (HIGH complexity; assess after VisiumHD adoption)

---

## Feature Prioritization Matrix

| Feature | User Value | Implementation Cost | Priority |
|---------|------------|---------------------|----------|
| Spatial plots: log1p_counts + log1p_genes + %mt, pre-filter | HIGH | LOW | P1 |
| Sample ordering by QC score in spatial section | HIGH | LOW | P1 |
| `qc_pass_spot` column in sct qc run | HIGH | MEDIUM | P1 |
| Pass/fail spatial overlay in pre-filter report | HIGH | LOW (given qc_pass_spot) | P1 |
| `sct concat` basic command with uns spatial merge | HIGH | MEDIUM | P1 |
| `sct concat` pipeline phase registration | HIGH | LOW | P1 |
| Celltype spatial overlay per sample (aggregate, categorical) | HIGH | MEDIUM | P1 |
| Consistent palette across samples for celltype spatial | HIGH | LOW | P1 |
| `subject_id` UMAP coloring in integration report | MEDIUM | LOW | P1 |
| Common color scale across samples (cohort vmin/vmax) | MEDIUM | MEDIUM | P2 |
| 4-panel per-sample spatial grid (2x2) | MEDIUM | MEDIUM | P2 |
| Split batch-vs-bio UMAP layout | MEDIUM | MEDIUM | P2 |
| `celltype_broad` spatial overlay | MEDIUM | LOW | P2 |
| Per-celltype spatial overlays (accordion, top-8) | MEDIUM | HIGH | P2 |
| `sct concat --from-dir` | LOW | LOW | P2 |
| Per-embedding UMAP colored by Leiden | LOW | MEDIUM | P3 |
| Spatial plots for IMC/Xenium | LOW | HIGH | P3 |
| Memory-efficient backed concat for VisiumHD | LOW | HIGH | P3 |

---

## Implementation Notes by Component

### `qc_pass_spot` in `sct qc run`
The per-sample `classify_samples` applies median-level thresholds. Per-spot pass/fail applies the same thresholds at spot level: `total_counts >= min_counts AND n_genes_by_counts >= min_genes AND pct_counts_mt <= max_pct_mt`. Before the `filter_cells` call, add `adata.obs['qc_pass_spot'] = True` then mark spots that would be removed as `False`. Store this in the h5ad so the pre-filter report can read it. The column name `qc_pass_spot` should be documented as a standard obs key.

### Spatial plot loop pattern
The pattern for pre-filter spatial plots follows what post-filter already does for `log1p_total_counts` (lines 302–334 in `report.py`). Extend to 4 metrics per sample. For the MVP, produce one figure per sample with multiple axes rather than one axis per metric per figure. Guard with `if str(lib) not in spatial_uns: continue` (already the pattern). Sort sample list by `classified['qc_pass']` ascending (failed first) then by `classified['total_counts_median']` ascending before the loop.

### `sct concat` uns merge correctness check
After `anndata.concat(..., uns_merge="unique")`, verify: `set(input_library_ids) == set(output.uns.get('spatial', {}).keys())`. If any library ID is missing from `uns['spatial']`, raise `SCToolsDataError` with a clear message listing which samples lost their spatial metadata. This is the single most important correctness invariant for downstream spatial plots.

### Celltype spatial section HTML size management
At 10 samples × 15 cell types = 150 plots at 150 DPI = approximately 1.5 MB per plot = 225 MB HTML. Browsers will hang or crash. MVP celltype spatial section uses **aggregate categorical overlay only** (one plot per sample showing all cell types). Per-celltype breakdown is a v1.x feature, implemented as an accordion with top-8-by-abundance limit.

---

## Sources

- Direct inspection: `sc_tools/qc/report.py` — confirmed which plots are generated per report type and which template slots are dead
- Direct inspection: `sc_tools/assets/pre_filter_qc_template.html`, `post_integration_qc_template.html`, `post_celltyping_qc_template.html` — confirmed template slot structure
- Direct inspection: `sc_tools/pl/spatial.py` — confirmed `plot_spatial_categorical` and `plot_spatial_continuous` exist with correct signatures
- Direct inspection: `sc_tools/qc/report_utils.py` — confirmed `auto_detect_embeddings`, `get_modality_terms`, `_CELL_MODALITIES` set
- Direct inspection: `sc_tools/cli/__init__.py` and `sc_tools/cli/` — confirmed no `concat` command exists
- Direct inspection: `.planning/PROJECT.md` — confirmed v2.0 milestone targets and Phase 6 subject metadata model status
- scanpy `sc.pl.spatial`: `ax=` parameter enables subplot integration; `vmin`/`vmax` enables common scale; `library_id` selects per-sample image from `uns['spatial']`
- AnnData `anndata.concat`: `uns_merge="unique"` preserves per-sample spatial metadata; `index_unique="-"` deduplicates obs names; `label` parameter injects batch key column

---

*Feature research for: sc_tools v2.0 Report Plots and Sample Concat milestone*
*Researched: 2026-03-27*
