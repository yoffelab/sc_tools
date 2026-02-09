# Project Architecture: Spatial Omics Analysis

This document outlines the directory structure and data flow for the computational oncology pipeline. All agentic operations must adhere to this roadmap to maintain project integrity and reproducibility.

## 1. Directory Overview

```text
.
├── Architecture.md      # System roadmap (this file)
├── skills.md            # Mandatory coding and statistical standards
├── Mission.md           # Active task tracking
├── Journal.md           # Technical and scientific decision log
├── environment.yml      # Conda environment configuration
├── data/                # Raw sequencing and imaging data (unprocessed)
├── metadata/            # Biological context (gene signatures, external obs)
├── results/             # Processed AnnData (.h5ad) and Seurat objects
├── scripts/             # Modular Python/R scripts for analysis
│   └── old code/        # Archived/Legacy script versions (Read-only)
├── figures/             # Output visualizations
│   ├── manuscript/      # Publication-ready plots (Strict statistical rules)
│   ├── spatial/         # Intermediate spatial plots of tissue slides
│   └── tls_spatial/     # Tertiary Lymphoid Structure specific analysis
├── outs/                # Platform-specific raw outputs (e.g., Visium/Xenium)
└── output/              # Intermediate tool outputs (e.g., Tangram batches)
```

---

## 2. Data Flow and Dependencies

The analysis follows a linear progression from data ingestion to high-level spatial statistics:

### Phase I: Ingestion & Preprocessing
* **Scripts:** `loupe2adata.py`, `annotation2mask_img.py`
* **Process:** Convert raw spatial formats to standardized `AnnData` or `SpatialData` containers. Clean up and store clinical variables such as solidity and architecture. Second script saves all images / masks within the anndata.
* **Output:** `results/results/adata.annotation.masks.h5ad`

### Phase II: QC, Integration & Generative Modeling
* **Scripts:** `preprocessing.py`
* **Process:** Preprocessing performs QC, cell / gene filtration, and data normalization. Mt, rb, and hb content is calculated, and hvg, svg are filtered at the end. Data is subsequently ran through **scVI** model on raw counts to learn latent representations and perform batch correction.
* **Output:** `results/scvi.h5ad`

### Phase III: Mapping & Deconvolution
* **Scripts:** `celltype_deconvolution_tangram.py`
* **Strategy:** Map single-cell reference data onto spatial coordinates using **Tangram** or **Cell2location**.
* **Output:** `results/adata.img.genescores.h5ad`

### Phase IV: Spatial Structures & Niche Analysis
* **Scripts:** `create_tls_anndata.py`, `tls_analysis.py`, `spot_deg.py`, `tumor_differences.py`, `macrophage_localization.py`, `process_colocalization.py`
* **Focus:** Identification of Tertiary Lymphoid Structures (TLS), EMT scoring, spatially variable gene (SVG) discovery, differential program analysis, macrophage colocalization, and process colocalization patterns.
* **Output:** `results/tls_clustered.h5ad`, `figures/manuscript/tumor_differences/`, `figures/manuscript/macrophage_localization/`, `figures/process_colocalization/`

### Phase V: Visualization & Figures
* **Scripts:** `manuscript_spatial_plots.py`, `score_gene_signatures.py`, `signature_heatmap.py`
* **Output:** `/figures/manuscript/`
* **Standard:** All plots must follow the statistical annotation rules (bars and asterisks) defined in `skills.md`.
* **Signature Heatmaps:** Comprehensive visualization of gene signature scores across spots with patient and solidity annotations, including both standard heatmaps and hierarchically clustered clustermaps.

---

## 3. Key Data Objects

| File Path | Description |
| :--- | :--- |
| `results/scvi.h5ad` | The primary integrated atlas with batch-corrected latent space. |
| `metadata/gene_signatures.json` | JSON dictionary of markers for Tumor, Stromal, and Immune programs. |
| `results/adata.annotation.masked.h5ad` | Final annotated object with specific cell-type masks. |
| `results/adata.img.genescores.h5ad` | AnnData with all gene signature scores (z-scored) for visualization and analysis. |
| `results/tls_clustered.h5ad` | Aggregation of TLS from all the tissue slides. |

---

## 4. Operational Rules for Agents

1. **Strict File Placement:** Never save plots in the root or `scripts/`. Always use the appropriate subfolder in `figures/`.
2. **Naming Convention:** New results should follow the format `adata.[step_description].h5ad` to maintain a clear audit trail.
3. **Statistical Rigor:** When generating boxplots or violin plots, apply Benjamini-Hochberg correction and draw significance asterisks as per `skills.md`. Numerical adjusted p-values must be included in the figure.
4. **Code Maintenance:** Do not modify `scripts/old code/`. If a legacy function is needed, refactor it into a new modular script.
5. **No Apostrophes:** Ensure all machine-generated documentation and emails avoid the use of apostrophes.

---

## 5. Development Environment
- **Environment:** Defined in `environment.yml` and `requirements.txt`.
- **Core Libraries:** `scvi-tools`, `spatialdata`, `squidpy`, `scanpy`, `statannotations`, `pinguoin`.