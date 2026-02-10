# sc_tools Pipeline Workflow

This document describes the **non-linear** analysis pipeline with human-in-the-loop phases. The workflow includes branching points and requires explicit input files (e.g., clinical metadata map) to bypass manual intervention.

**Rendering:** The Mermaid diagram below renders in GitHub, GitLab, and many Markdown viewers. To export as PNG/SVG, use [Mermaid CLI](https://github.com/mermaid-js/mermaid-cli): `mmdc -i WORKFLOW.md -o workflow_diagram.png`.

---

## Workflow Diagram

```mermaid
flowchart TB
    subgraph P1["Phase 1: Data Ingestion & QC"]
        direction TB
        A1[Platform-specific ingestion] --> A2[Raw unnormalized AnnData]
        A2 --> A3["QC metrics (scanpy.pp.calculate_qc_metrics)"]
        A3 --> A4["QC report: figures/QC/raw/"]
        A1 --> |"Visium/HD: fastq + H&E + Cytassist → spaceranger → cloupe"| A2
        A1 --> |"IMC: mcd/txt → segmentation → h5ad"| A2
        A1 --> |"Xenium: assume preprocessed"| A2
    end

    subgraph P2["Phase 2: Metadata Attachment"]
        direction TB
        B1{Clinical map file provided?} --> |"No: HUMAN-IN-LOOP"| B2[Prepare sample→metadata CSV/xlsx]
        B2 --> B1
        B1 --> |"Yes: $(PROJECT)/metadata/sample_metadata.csv"| B3[Join clinical metadata to adata.obs]
        B3 --> B4[All clinical metadata attached]
    end

    subgraph P3["Phase 3: Preprocessing"]
        direction TB
        C1[Backup adata.raw] --> C2[Filter cells/genes/samples failing QC]
        C2 --> C3[Normalize, batch correct, cluster]
        C3 --> C4[Automated cell typing]
        C4 --> C5["QC report: figures/QC/post/"]
    end

    subgraph P35["Phase 3.5: Demographics (branch)"]
        D1[Cohort stats: piechart, violin, bar, heatmap, etc.]
        D2[Figure 1 for manuscript]
    end

    subgraph P4["Phase 4: Manual Cell Typing"]
        direction TB
        E1[Extract cluster_id from adata.obs] --> E2[Provide JSON: cluster_id→celltype]
        E2 --> E3{Mapping satisfactory?}
        E3 --> |"No: HUMAN-IN-LOOP"| E2
        E3 --> |"Yes"| E4[Apply celltype + celltype_broad]
        E4 --> E5[Matrixplot, UMAP, signature genes]
    end

    subgraph P5["Phase 5: Downstream Biology"]
        F1[Gene scoring, deconvolution]
        F2[Spatial/process analysis]
        F3[Colocalization, Moran's I, neighborhood enrichment]
    end

    subgraph P67["Phase 6–7: Meta Analysis (optional)"]
        G1[Phase 6: Aggregate to ROI + patient level]
        G2[Phase 7: Downstream on aggregated data]
    end

    P1 --> P2
    P2 --> P3
    P3 --> P35
    P3 --> P4
    P4 --> P5
    P3 --> |"Skip Phase 4 if auto-typing sufficient"| P5
    P5 --> P67

    P2 -.-> |"Preprocessed projects may START here"| P3
    P3 -.-> |"Already clustered projects may START here"| P4
```

---

## Phase Summary

| Phase | Name | Human-in-Loop? | Required Input | Output |
|-------|------|----------------|----------------|--------|
| **1** | Data Ingestion & QC | No | Platform raw data | Raw AnnData, `$(PROJECT)/figures/QC/raw/` |
| **2** | Metadata Attachment | Yes (unless map provided) | `$(PROJECT)/metadata/sample_metadata.csv` or `.xlsx` | AnnData with clinical metadata |
| **3** | Preprocessing | No | — | Filtered, normalized, clustered AnnData; `$(PROJECT)/figures/QC/post/` |
| **3.5** | Demographics | Project-specific | — | Figure 1, cohort stats |
| **4** | Manual Cell Typing | Yes (iterative) | JSON: `cluster_id→celltype` | Phenotyped AnnData |
| **5** | Downstream Biology | No | — | Gene scores, spatial analysis, figures |
| **6–7** | Meta Analysis | No | — | ROI/patient aggregated results |

---

## Key File Descriptions

All paths are project-specific under `projects/<platform>/<project_name>/`.

| File / Path | Description |
|-------------|-------------|
| `metadata/sample_metadata.csv` or `.xlsx` | Sample→clinical metadata map. Enables Phase 2 without human-in-loop. Columns: `sample` (matches `adata.obs['sample']`), plus any clinical columns. |
| `metadata/celltype_map.json` | Cluster→celltype mapping for Phase 4. Format: `{cluster_id: {celltype_name: "...", total_obs_count: N}}`. cluster_id type must match `adata.obs` (string/int). |
| `adata.obs['sample']` | Annotates original sample origin. Required after Phase 1. |
| `adata.obs['raw_data_dir']` | Backup path for original data location. Required after Phase 1. |
| `adata.obsm['spatial']` | Spatial coordinates for spatial platforms. |
| `figures/QC/raw/` | Pre-normalization QC reports (2x2 histogram grid, spatial multipage). |
| `figures/QC/post/` | Post-normalization QC reports. |
