# Plan: sc_tools Toolkit and Pipeline

**Scope:** Repository-level and generalizable goals. Project-specific plans live in each project's `Plan.md` under `projects/<data_type>/<project_name>/Plan.md`.

**Last Updated:** 2026-03-12

---

## Pipeline Implementation Status

Pipeline phases and checkpoint specs: see [Architecture.md](Architecture.md) §2.

### ingest_raw / ingest_load

- [x] Batch manifest system (`sc_tools.ingest.config`)
- [x] Space Ranger command builder (`sc_tools.ingest.spaceranger`)
- [x] Xenium Ranger command builder (`sc_tools.ingest.xenium`)
- [x] IMC pipeline command builder (`sc_tools.ingest.imc`)
- [x] Snakemake `ingest_raw` rules in template
- [x] Modality loaders: Visium, Visium HD, Visium HD Cell, Xenium, IMC
- [ ] CosMx loader (deprioritized)
- [x] Visium HD Cell modality (SpaceRanger 4 cell segmentation)
- [x] Concatenation: `concat_samples()`
- [ ] IMC end-to-end Snakemake workflow for lymph_dlbcl and ggo-imc
- [x] IMC image loading (`load_imc_sample(load_images=True)`)
- [x] H&E image loading (any modality)
- [x] IMC spatial plotting (`plot_imc_composite`, `plot_imc_channel`)
- [x] Generic `scripts/ingest.py` checkpoint script
- [ ] `ingest_load` Snakemake rule with sentinel
- [ ] SpatialData (optional, Visium HD / Xenium)

### Checkpoint Validation

- [x] `sc_tools.validate` — validate_p1 through validate_p4, auto-fix
- [x] CLI: `scripts/validate_checkpoint.py`
- [x] Snakemake validation sentinels

### qc_filter

- [x] Load `ingest_load` checkpoints
- [x] Per-sample QC: `filter_spots()`, modality-aware thresholds
- [x] Concatenation across samples
- [x] QC metrics, plots, HTML reports (4 date-versioned reports)
- [x] Sample-level QC (`classify_samples`, `apply_qc_filter`)
- [ ] MA plots (future)

### metadata_attach

- [ ] Bypass file support (`sample_metadata.csv` / `.xlsx`)
- [ ] Human-guided join when no file provided

### preprocess

- [x] `sc_tools.pp` module (modality-aware, GPU auto-detection)
- [x] Normalization, integration (scVI, Harmony, CytoVI), clustering
- [x] Recipes: `preprocess(modality=...)`
- [x] Integration benchmark workflow
- [x] 4 QC HTML reports
- [ ] resolVI integration
- [ ] Post-normalization QC report

### demographics

- [ ] sc_tools helpers: piechart, histogram, violin, bar, heatmap
- [ ] Figure 1 for cohort description

### scoring

- [x] Gene signature scoring (scanpy, ucell, ssgsea)
- [x] Gene set loaders (Hallmark bundled, MSigDB, GMT)
- [x] Group-level enrichment (ORA, GSEA pseudobulk)
- [x] Automated cell typing (`sc_tools.tl.celltype`)
- [x] Cell-type deconvolution (cell2location, tangram, destvi)

### celltype_manual

- [ ] sc_tools workflow: extract clusters, JSON template
- [ ] Iterative refinement with matrixplot, UMAP, scatter
- [x] Post-celltyping QC report

### biology / meta_analysis

- [ ] Spatial/process analysis, colocalization, neighborhood enrichment
- [ ] ROI/patient aggregation

---

## Testing

- [ ] ggo_visium project tests
- [ ] sc_tools package test expansion
- Implementation order: project tests → sc_tools tests → implement

---

## CI/CD Roadmap

| Step | Status |
|------|--------|
| Linting (Ruff) | Done |
| Snakemake | Done |
| Sphinx docs | Done |
| PyPI deployment | Done |
| GitHub Actions | Done |
| Storage + Registry + MCP | Done |
| Registry + Phase DAG | Done |
| BioData v1.0 | Done |
| BioData v1.1 | Done |

---

## To Do (later)

- [ ] Organize production scripts by phase within project
- [ ] Update imports to use `sc_tools.*` everywhere
- [ ] CosMx loader for 1k/6k/full_library panels
- [ ] Documentation: migration guide, docstrings, API docs

---

## Reference

- **README.md** — Pipeline workflow diagram, phase summary
- **Architecture.md** — Directory layout, phase details, project paths, testing
- **skills.md** — Coding and statistical standards
- **Per-project:** `projects/<data_type>/<project_name>/Plan.md` and `Journal.md`
