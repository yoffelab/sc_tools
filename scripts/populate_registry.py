#!/usr/bin/env python
"""Populate the sc_tools registry with all known local and public datasets.

Run from repo root:
    python scripts/populate_registry.py
"""

from __future__ import annotations

import os
from pathlib import Path

BASE = Path("/Users/junbumkim/Documents")

# ── Shared registry connection ────────────────────────────────────────────────
from sc_tools.registry import Registry

r = Registry()


def add_project_safe(name: str, **kwargs) -> None:
    """Add a project (idempotent)."""
    try:
        r.add_project(name, **kwargs)
        print(f"  [+] project: {name}")
    except Exception as e:
        print(f"  [~] project exists or error ({name}): {e}")


def reg(
    project: str,
    phase: str,
    path: str | Path,
    *,
    file_role: str = "primary",
    sample_id: str | None = None,
    entry_phase: bool = False,
    status: str = "ready",
    n_obs: int | None = None,
    n_vars: int | None = None,
    n_samples: int | None = None,
    notes: str | None = None,
    fmt: str = "h5ad",
) -> None:
    """Register a dataset and upsert its phase row.  Skips missing local files."""
    uri = str(path)
    if uri.startswith("/") and not Path(uri).exists():
        print(f"  [!] skip (not found): {uri}")
        return
    try:
        ds_id = r.register_dataset(
            project,
            phase,
            uri,
            sample_id=sample_id,
            fmt=fmt,
            file_role=file_role,
            status=status,
        )
        r.upsert_phase(
            project,
            phase,
            status="complete",
            entry_phase=entry_phase,
            primary_dataset_id=ds_id if file_role in ("primary", "entry_point") else None,
            n_obs=n_obs,
            n_vars=n_vars,
            n_samples=n_samples,
            notes=notes,
        )
        print(f"  [+] {project} / {phase}: {Path(uri).name}")
    except Exception as e:
        print(f"  [~] already exists or error ({project}/{phase}): {e}")


def reg_pending(project: str, phase: str, url: str, notes: str | None = None) -> None:
    """Register a public-source dataset (URI = download URL, status=pending)."""
    try:
        ds_id = r.register_dataset(
            project,
            phase,
            url,
            status="pending",
            file_role="entry_point",
        )
        r.upsert_phase(
            project,
            phase,
            status="not_started",
            entry_phase=True,
            primary_dataset_id=ds_id,
            notes=notes,
        )
        print(f"  [+] {project} / {phase}: (pending)")
    except Exception as e:
        print(f"  [~] {project}/{phase}: {e}")


# ==============================================================================
# A. LOCAL PROJECTS
# ==============================================================================
print("\n" + "=" * 70)
print("A. LOCAL PROJECTS")
print("=" * 70)

# ── imc/aapc ─────────────────────────────────────────────────────────────────
print("\n--- imc/aapc ---")
add_project_safe(
    "aapc",
    platform="imc",
    data_type="imc",
    domain="spatial_proteomics",
    imaging_modality="mass_spec_imaging",
    project_type="internal",
    visibility="private",
)
P = BASE / "sc_tools/projects/imc/aapc"
for panel in ("panel1", "panel2"):
    reg("aapc", "qc_filter",        P / f"results/{panel}.raw.h5ad",         sample_id=panel, notes=f"AAPC {panel} raw")
    reg("aapc", "preprocess",       P / f"results/{panel}.harmony.h5ad",     sample_id=panel)
    reg("aapc", "celltype_manual",  P / f"results/{panel}.clin_celltyped.h5ad", sample_id=panel)
    reg("aapc", "biology",          P / f"results/{panel}.roi.celltype.h5ad", sample_id=panel, fmt="h5ad")

# ── imc/bcg_bladder ──────────────────────────────────────────────────────────
print("\n--- imc/bcg_bladder ---")
add_project_safe(
    "bcg_bladder",
    platform="imc",
    data_type="imc",
    domain="spatial_proteomics",
    imaging_modality="mass_spec_imaging",
    project_type="internal",
    visibility="private",
)
P = BASE / "sc_tools/projects/imc/bcg_bladder"
reg("bcg_bladder", "qc_filter",  P / "results/quant.raw.h5ad")
reg("bcg_bladder", "preprocess", P / "results/quant.harmony.leiden.h5ad")

# ── imc/brca-unet ────────────────────────────────────────────────────────────
print("\n--- imc/brca-unet ---")
add_project_safe(
    "brca_unet",
    platform="imc",
    data_type="imc",
    domain="spatial_proteomics",
    imaging_modality="mass_spec_imaging",
    project_type="internal",
    visibility="private",
)
P = BASE / "sc_tools/projects/imc/brca-unet"
reg("brca_unet", "preprocess", P / "adata.log.scale.h5ad", notes="Pre-integrated log-scale; non-standard path")

# ── imc/clonevo ──────────────────────────────────────────────────────────────
print("\n--- imc/clonevo ---")
add_project_safe(
    "clonevo",
    platform="imc",
    data_type="imc",
    domain="spatial_proteomics",
    imaging_modality="mass_spec_imaging",
    project_type="internal",
    visibility="private",
)
P = BASE / "sc_tools/projects/imc/clonevo"
reg("clonevo", "qc_filter",       P / "results/raw.h5ad")
reg("clonevo", "preprocess",      P / "results/quant.harmony.h5ad")
reg("clonevo", "celltype_manual", P / "results/phenotyped.h5ad")
reg("clonevo", "celltype_manual", P / "results/phenotyped.labeled.h5ad",
    file_role="supplementary", notes="Labeled variant")
for subtype in ("fibs", "immune", "lymphocytes", "tumor", "stromal"):
    reg("clonevo", "biology", P / f"results/{subtype}.h5ad", sample_id=subtype,
        file_role="supplementary")

# ── imc/colon_liver_meta ─────────────────────────────────────────────────────
print("\n--- imc/colon_liver_meta ---")
add_project_safe(
    "colon_liver_meta",
    platform="imc",
    data_type="imc",
    domain="spatial_proteomics",
    imaging_modality="mass_spec_imaging",
    project_type="internal",
    visibility="private",
)
# No h5ad files found; project stub only

# ── imc/ggo-imc ──────────────────────────────────────────────────────────────
print("\n--- imc/ggo-imc ---")
# Already exists as ggo_human; check
P = BASE / "sc_tools/projects/imc/ggo-imc"
for panel in ("panel_g", "panel_h"):
    reg("ggo_human", "preprocess",      P / f"results/{panel}.harmony.h5ad", sample_id=panel)
    reg("ggo_human", "celltype_manual", P / f"results/{panel}.celltyped.h5ad", sample_id=panel)
    reg("ggo_human", "biology",         P / f"results/{panel}.roi.celltype.h5ad", sample_id=panel)

# ── imc/healthy_lung ─────────────────────────────────────────────────────────
print("\n--- imc/healthy_lung ---")
add_project_safe(
    "healthy_lung",
    platform="imc",
    data_type="imc",
    domain="spatial_proteomics",
    imaging_modality="mass_spec_imaging",
    project_type="internal",
    visibility="private",
)
P = BASE / "sc_tools/projects/imc/healthy_lung"
reg("healthy_lung", "qc_filter",       P / "IMC_healthy_lung_preprocessed_no_COPD.h5ad",
    notes="Pre-processed; non-standard path")
reg("healthy_lung", "preprocess",      P / "IMC_healthy_lung_batch_corrected.h5ad",
    notes="Batch corrected; non-standard path")
reg("healthy_lung", "celltype_manual", P / "gatdu_labelled.h5ad",
    notes="Manually labeled; non-standard path")

# ── imc/imc_utuc ─────────────────────────────────────────────────────────────
print("\n--- imc/imc_utuc ---")
add_project_safe(
    "imc_utuc",
    platform="imc",
    data_type="imc",
    domain="spatial_proteomics",
    imaging_modality="mass_spec_imaging",
    project_type="internal",
    visibility="private",
)
P = BASE / "sc_tools/projects/imc/imc_utuc"
reg("imc_utuc", "preprocess", P / "UTUC_adata.h5ad", notes="Pre-processed; non-standard path")

# ── imc/pancreas_liver_meta ───────────────────────────────────────────────────
print("\n--- imc/pancreas_liver_meta ---")
add_project_safe(
    "pancreas_liver_meta",
    platform="imc",
    data_type="imc",
    domain="spatial_proteomics",
    imaging_modality="mass_spec_imaging",
    project_type="internal",
    visibility="private",
)
P = BASE / "sc_tools/projects/imc/pancreas_liver_meta"
reg("pancreas_liver_meta", "qc_filter",       P / "results/quant.raw.h5ad")
reg("pancreas_liver_meta", "preprocess",      P / "results/quant.harmony.leiden.h5ad")
reg("pancreas_liver_meta", "celltype_manual", P / "results/phenotyped.h5ad")
reg("pancreas_liver_meta", "biology",         P / "results/patient_celltype.h5ad", sample_id="patient_agg")
reg("pancreas_liver_meta", "biology",         P / "results/roi_celltype.h5ad",     sample_id="roi_agg",
    file_role="supplementary")

# ── visium_hd/10xCRC ─────────────────────────────────────────────────────────
print("\n--- visium_hd/10xCRC ---")
add_project_safe(
    "10xCRC",
    platform="visium_hd",
    data_type="visium_hd",
    domain="spatial_transcriptomics",
    imaging_modality="sequencing_based",
    project_type="internal",
    visibility="private",
)
P = BASE / "sc_tools/projects/visium_hd/10xCRC"
reg("10xCRC", "qc_filter",  P / "output/raw_data.h5ad",          notes="Non-standard output/ path")
reg("10xCRC", "qc_filter",  P / "output/adata_filtered.h5ad",    file_role="supplementary")
reg("10xCRC", "preprocess", P / "output/adata_corrected.h5ad")
reg("10xCRC", "preprocess", P / "output/adata_corrected_scvi.h5ad",
    file_role="supplementary", notes="scVI integration variant")

# ── cosmx_1k/lymph_dlbcl ─────────────────────────────────────────────────────
print("\n--- cosmx_1k/lymph_dlbcl ---")
add_project_safe(
    "cosmx_lymph_dlbcl",
    platform="cosmx",
    data_type="cosmx_1k",
    domain="spatial_transcriptomics",
    imaging_modality="probe_based",
    project_type="internal",
    visibility="private",
)
# No h5ad files found; project stub only

# ── visium_hd_cell/robin ──────────────────────────────────────────────────────
print("\n--- visium_hd_cell/robin ---")
add_project_safe(
    "robin_cell",
    platform="visium_hd_cell",
    data_type="visium_hd_cell",
    domain="spatial_transcriptomics",
    imaging_modality="sequencing_based",
    project_type="internal",
    visibility="private",
)
# Empty results; project stub only


# ==============================================================================
# B. PUBLIC / EXTERNAL DATASETS
# ==============================================================================
print("\n" + "=" * 70)
print("B. PUBLIC / EXTERNAL DATASETS")
print("=" * 70)


def pub(name: str, platform: str, data_type: str, domain: str, imaging_modality: str,
        url: str, notes: str | None = None) -> None:
    """Add an external public dataset as a project stub with pending entry_point."""
    add_project_safe(
        name,
        platform=platform,
        data_type=data_type,
        domain=domain,
        imaging_modality=imaging_modality,
        project_type="external",
        visibility="public",
    )
    reg_pending(name, "ingest_raw", url, notes=notes)


# ── 10x Genomics Visium Classic ───────────────────────────────────────────────
print("\n--- 10x Visium (classic) ---")
VISIUM_BASE = "https://www.10xgenomics.com/datasets"

pub("10x_visium_human_breast_cancer",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/visium-hd-cytassist-gene-expression-libraries-of-post-xenium-human-breast-cancer",
    "Human breast cancer; 10x Genomics public dataset")

pub("10x_visium_human_heart",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/spatial-gene-expression-human-heart",
    "Human heart spatial gene expression")

pub("10x_visium_human_brain_coronal",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/spatial-gene-expression-human-cerebral-cortex",
    "Human cerebral cortex")

pub("10x_visium_human_ovarian_cancer",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/spatial-gene-expression-v1-human-ovarian-cancer",
    "Human ovarian cancer")

pub("10x_visium_mouse_brain_coronal",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/spatial-gene-expression-mouse-brain-coronal-section-1",
    "Mouse brain coronal section 1")

pub("10x_visium_mouse_brain_sagittal_ant",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/spatial-gene-expression-mouse-brain-sagittal-anterior",
    "Mouse brain sagittal anterior")

pub("10x_visium_mouse_brain_sagittal_post",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/spatial-gene-expression-mouse-brain-sagittal-posterior",
    "Mouse brain sagittal posterior")

pub("10x_visium_mouse_kidney",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/spatial-gene-expression-mouse-kidney",
    "Mouse kidney")

pub("10x_visium_human_lymph_node",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/spatial-gene-expression-human-lymph-node",
    "Human lymph node")

pub("10x_visium_human_skin",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/spatial-gene-expression-human-skin",
    "Human skin")

pub("10x_visium_human_colon_cancer",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/spatial-gene-expression-human-colon-cancer-1",
    "Human colon cancer")

pub("10x_visium_human_glioblastoma",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/spatial-gene-expression-human-glioblastoma",
    "Human glioblastoma")

pub("10x_visium_human_prostate_cancer",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/spatial-gene-expression-human-prostate-cancer",
    "Human prostate cancer")

pub("10x_visium_cytassist_human_lung_cancer",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/visium-cytassist-gene-expression-human-lung-cancer",
    "CytAssist human lung cancer")

pub("10x_visium_cytassist_human_melanoma",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/visium-cytassist-gene-expression-human-melanoma",
    "CytAssist human melanoma")

pub("10x_visium_cytassist_human_crc",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/visium-cytassist-gene-expression-human-colorectal-cancer",
    "CytAssist human colorectal cancer")

pub("10x_visium_cytassist_mouse_brain",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/visium-cytassist-gene-expression-mouse-brain",
    "CytAssist mouse brain")

pub("10x_visium_cytassist_human_cervical_cancer",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/visium-cytassist-gene-expression-human-cervical-cancer",
    "CytAssist human cervical cancer")

pub("10x_visium_cytassist_human_kidney",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/visium-cytassist-gene-expression-human-kidney",
    "CytAssist human kidney")

pub("10x_visium_cytassist_human_tonsil",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/visium-cytassist-gene-expression-human-tonsil",
    "CytAssist human tonsil")

pub("10x_visium_ffpe_human_breast_cancer",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/spatial-gene-expression-v1-human-breast-cancer-ffpe",
    "FFPE human breast cancer; archival tissue")

pub("10x_visium_ffpe_human_prostate",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/spatial-gene-expression-v1-human-prostate-ffpe",
    "FFPE human prostate tissue")

pub("10x_visium_ffpe_human_cerebellum",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/spatial-gene-expression-human-cerebellum-ffpe",
    "FFPE human cerebellum")

# ── 10x Visium HD ─────────────────────────────────────────────────────────────
print("\n--- 10x Visium HD ---")

pub("10x_visium_hd_human_crc",
    "visium_hd", "visium_hd", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/visium-hd-cytassist-gene-expression-libraries-of-human-crc",
    "Visium HD human colorectal cancer; 2um bin resolution")

pub("10x_visium_hd_human_breast_cancer",
    "visium_hd", "visium_hd", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/visium-hd-cytassist-gene-expression-libraries-of-human-breast-cancer",
    "Visium HD human breast cancer")

pub("10x_visium_hd_mouse_brain",
    "visium_hd", "visium_hd", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain",
    "Visium HD mouse brain")

pub("10x_visium_hd_mouse_intestine",
    "visium_hd", "visium_hd", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/visium-hd-cytassist-gene-expression-libraries-of-mouse-intestine",
    "Visium HD mouse intestine")

pub("10x_visium_hd_human_prostate_cancer",
    "visium_hd", "visium_hd", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/visium-hd-cytassist-gene-expression-libraries-of-human-prostate-cancer",
    "Visium HD human prostate cancer")

pub("10x_visium_hd_human_pancreatic_cancer",
    "visium_hd", "visium_hd", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/visium-hd-cytassist-gene-expression-libraries-of-human-pancreatic-cancer",
    "Visium HD human pancreatic cancer; paired with Xenium")

pub("10x_visium_hd_post_xenium_human_bc",
    "visium_hd", "visium_hd", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/visium-hd-cytassist-gene-expression-libraries-of-post-xenium-human-breast-cancer",
    "Visium HD post-Xenium human breast cancer; multi-modal pairing dataset")

pub("10x_visium_hd_human_lung_cancer",
    "visium_hd", "visium_hd", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/visium-hd-cytassist-gene-expression-libraries-of-human-lung-cancer",
    "Visium HD human lung cancer")

pub("10x_visium_hd_human_crc_spacerangerhd4",
    "visium_hd", "visium_hd", "spatial_transcriptomics", "sequencing_based",
    f"{VISIUM_BASE}/visium-hd-space-ranger-4-human-crc",
    "Visium HD SpaceRanger 4.0 CRC; cell segmentation output")

# ── 10x Xenium ────────────────────────────────────────────────────────────────
print("\n--- 10x Xenium ---")
XENIUM_BASE = "https://www.10xgenomics.com/datasets"

pub("10x_xenium_human_breast_cancer",
    "xenium", "xenium", "spatial_transcriptomics", "probe_based",
    f"{XENIUM_BASE}/xenium-human-breast-cancer-gene-expression",
    "Xenium human breast cancer; 2 replicates; 280 gene panel")

pub("10x_xenium_human_lung_cancer",
    "xenium", "xenium", "spatial_transcriptomics", "probe_based",
    f"{XENIUM_BASE}/xenium-human-lung-cancer-gene-expression",
    "Xenium human NSCLC; 313 gene panel")

pub("10x_xenium_human_brain",
    "xenium", "xenium", "spatial_transcriptomics", "probe_based",
    f"{XENIUM_BASE}/xenium-human-brain-gene-expression",
    "Xenium human brain (frontal cortex + hippocampus); 392 gene panel")

pub("10x_xenium_human_skin",
    "xenium", "xenium", "spatial_transcriptomics", "probe_based",
    f"{XENIUM_BASE}/xenium-human-skin-gene-expression",
    "Xenium human skin; 377 gene panel")

pub("10x_xenium_human_colon_cancer",
    "xenium", "xenium", "spatial_transcriptomics", "probe_based",
    f"{XENIUM_BASE}/xenium-human-colon-cancer-gene-expression",
    "Xenium human colon cancer")

pub("10x_xenium_human_pancreatic_cancer",
    "xenium", "xenium", "spatial_transcriptomics", "probe_based",
    f"{XENIUM_BASE}/xenium-human-pancreatic-cancer-gene-expression",
    "Xenium human pancreatic cancer; paired with Visium HD")

pub("10x_xenium_human_lymph_node",
    "xenium", "xenium", "spatial_transcriptomics", "probe_based",
    f"{XENIUM_BASE}/xenium-human-lymph-node-gene-expression",
    "Xenium human lymph node")

pub("10x_xenium_human_prostate",
    "xenium", "xenium", "spatial_transcriptomics", "probe_based",
    f"{XENIUM_BASE}/xenium-human-prostate-gene-expression",
    "Xenium human prostate")

pub("10x_xenium_mouse_brain",
    "xenium", "xenium", "spatial_transcriptomics", "probe_based",
    f"{XENIUM_BASE}/xenium-mouse-brain-gene-expression",
    "Xenium mouse brain (adult coronal; 1000+ panel)")

pub("10x_xenium_mouse_liver",
    "xenium", "xenium", "spatial_transcriptomics", "probe_based",
    f"{XENIUM_BASE}/xenium-mouse-liver-gene-expression",
    "Xenium mouse liver")

pub("10x_xenium_human_crc_5k",
    "xenium", "xenium", "spatial_transcriptomics", "probe_based",
    f"{XENIUM_BASE}/xenium-human-colorectal-cancer-gene-expression-5000-panel",
    "Xenium human CRC; 5000-gene panel")

pub("10x_xenium_human_ovarian_cancer",
    "xenium", "xenium", "spatial_transcriptomics", "probe_based",
    f"{XENIUM_BASE}/xenium-human-ovarian-cancer-gene-expression",
    "Xenium human ovarian cancer")

pub("10x_xenium_human_melanoma",
    "xenium", "xenium", "spatial_transcriptomics", "probe_based",
    f"{XENIUM_BASE}/xenium-human-melanoma-gene-expression",
    "Xenium human melanoma")

# ── 10x Chromium (scRNA / multiome) ──────────────────────────────────────────
print("\n--- 10x Chromium (scRNA reference datasets) ---")
CHROMIUM_BASE = "https://www.10xgenomics.com/datasets"

pub("10x_pbmc_3k",
    "10x_chromium", "scrna", "single_cell", "sequencing_based",
    f"{CHROMIUM_BASE}/pbmc-3k-filtered-feature-bc-matrix",
    "PBMC 3k; canonical scRNA benchmark dataset")

pub("10x_pbmc_10k_v3",
    "10x_chromium", "scrna", "single_cell", "sequencing_based",
    f"{CHROMIUM_BASE}/10k-human-pbmcs-3-v3-1-chromium-x",
    "PBMC 10k; Chromium X v3.1")

pub("10x_pbmc_multiome",
    "10x_chromium", "multiome", "single_cell", "sequencing_based",
    f"{CHROMIUM_BASE}/pbmc-granulocyte-sorted-10k-filtered-feature-bc-matrix",
    "PBMC 10k multiome (RNA+ATAC)")

pub("10x_human_breast_cancer_multiome",
    "10x_chromium", "multiome", "single_cell", "sequencing_based",
    f"{CHROMIUM_BASE}/fresh-frozen-human-breast-cancer-1",
    "Human breast cancer multiome (RNA+ATAC)")

pub("10x_human_brain_3k_multiome",
    "10x_chromium", "multiome", "single_cell", "sequencing_based",
    f"{CHROMIUM_BASE}/human-brain-3k-multi-ome",
    "Human brain 3k multiome (RNA+ATAC)")

# ── IMC public datasets (imcdatasets Python package) ─────────────────────────
print("\n--- IMC public datasets (imcdatasets) ---")
IMC_DATASETS = "https://zenodo.org/record"

pub("imc_damond2019_pancreas",
    "imc", "imc", "spatial_proteomics", "mass_spec_imaging",
    "https://zenodo.org/record/3951011",
    "Damond et al. 2019 Cell Metabolism; T1D pancreas; 35 markers; available via imcdatasets pkg")

pub("imc_jackson2020_breast_cancer",
    "imc", "imc", "spatial_proteomics", "mass_spec_imaging",
    "https://zenodo.org/record/7575859",
    "Jackson & Fischer et al. 2020 Nature; breast cancer TME; 38 markers; imcdatasets pkg")

pub("imc_keren2018_tnbc",
    "imc", "imc", "spatial_proteomics", "mass_spec_imaging",
    "https://www.angelolab.com/mibi-data",
    "Keren et al. 2018 Cell; TNBC; MIBI-TOF (similar to IMC); 36 markers")

pub("imc_schapiro2022_lung_crc",
    "imc", "imc", "spatial_proteomics", "mass_spec_imaging",
    "https://zenodo.org/record/5949116",
    "Schapiro et al. 2022 Nature Methods; lung + CRC; 9 channels; tutorial dataset")

pub("imc_hoch2022_melanoma",
    "imc", "imc", "spatial_proteomics", "mass_spec_imaging",
    "https://zenodo.org/record/7637156",
    "Hoch et al. 2022 Science Immunology; melanoma; 35 markers")

pub("imc_crc_atlas",
    "imc", "imc", "spatial_proteomics", "mass_spec_imaging",
    "https://zenodo.org/record/6024273",
    "Hartmann et al. 2020 biorXiv; CRC atlas; 35 markers; imcdatasets pkg (IMMUcan_2022_CancerExample)")

pub("imc_basel_zurich_breast",
    "imc", "imc", "spatial_proteomics", "mass_spec_imaging",
    "https://zenodo.org/record/3518284",
    "Ali et al. 2020 Nature Cancer; breast cancer; Basel+Zurich cohorts; 38 markers")

pub("imc_wagner2019_breast_cancer",
    "imc", "imc", "spatial_proteomics", "mass_spec_imaging",
    "https://data.mendeley.com/datasets/yz3s8fgxmj",
    "Wagner et al. 2019 Cell; breast cancer; 73-marker CyTOF + MIBI")

pub("imc_covid19_lung",
    "imc", "imc", "spatial_proteomics", "mass_spec_imaging",
    "https://zenodo.org/record/5018260",
    "Rendeiro et al. 2021 Nature; COVID-19 lung autopsy; 37 markers")

pub("imc_human_spleen",
    "imc", "imc", "spatial_proteomics", "mass_spec_imaging",
    "https://zenodo.org/record/5949116",
    "Catena et al. 2023; human spleen; 32 markers; imcdatasets example data")

# ── MERFISH / Allen Brain / AIBS ──────────────────────────────────────────────
print("\n--- MERFISH / Allen Brain ---")
pub("merfish_mouse_brain_abc_atlas",
    "merfish", "merfish", "spatial_transcriptomics", "probe_based",
    "https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas",
    "Allen Brain Cell Atlas; mouse whole-brain MERFISH; ~4M cells; 500 gene panel")

pub("merfish_human_brain_abc_atlas",
    "merfish", "merfish", "spatial_transcriptomics", "probe_based",
    "https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas",
    "Allen Brain Cell Atlas; human multi-region MERFISH")

pub("merfish_mouse_mop",
    "merfish", "merfish", "spatial_transcriptomics", "probe_based",
    "https://zenodo.org/record/8167488",
    "AIBS mouse primary motor cortex MERFISH; 254 genes; published in Science 2021")

pub("merfish_mouse_hypothalamus",
    "merfish", "merfish", "spatial_transcriptomics", "probe_based",
    "https://datadryad.org/stash/dataset/doi:10.5061/dryad.8t8s248",
    "Moffitt et al. 2018 Science; mouse hypothalamic preoptic region; 155 genes")

pub("merfish_human_ovarian_cancer",
    "merfish", "merfish", "spatial_transcriptomics", "probe_based",
    "https://zenodo.org/record/7957748",
    "MERFISH human ovarian cancer TME; 300 gene panel")

pub("merfish_mouse_ileum",
    "merfish", "merfish", "spatial_transcriptomics", "probe_based",
    "https://zenodo.org/record/7319872",
    "MERFISH mouse ileum; 241 genes; Gut Cell Atlas subset")

# ── Vizgen MERSCOPE ────────────────────────────────────────────────────────────
print("\n--- Vizgen MERSCOPE ---")
pub("vizgen_merscope_human_liver_cancer",
    "merscope", "merscope", "spatial_transcriptomics", "probe_based",
    "https://console.cloud.google.com/storage/browser/vz-liver-showcase",
    "Vizgen liver cancer showcase; 500 gene panel; ~2 million cells")

pub("vizgen_merscope_human_crc",
    "merscope", "merscope", "spatial_transcriptomics", "probe_based",
    "https://info.vizgen.com/merscope-data-release-human-crc",
    "Vizgen MERSCOPE human colorectal cancer data release; 500 genes")

pub("vizgen_merscope_mouse_brain",
    "merscope", "merscope", "spatial_transcriptomics", "probe_based",
    "https://console.cloud.google.com/storage/browser/vz-brain-showcase",
    "Vizgen mouse brain showcase; 483 genes; coronal + sagittal sections")

pub("vizgen_merscope_human_breast_cancer",
    "merscope", "merscope", "spatial_transcriptomics", "probe_based",
    "https://info.vizgen.com/merscope-data-release-human-breast-cancer",
    "Vizgen MERSCOPE human breast cancer data release")

pub("vizgen_merscope_human_kidney",
    "merscope", "merscope", "spatial_transcriptomics", "probe_based",
    "https://console.cloud.google.com/storage/browser/vz-kidney-showcase",
    "Vizgen kidney showcase dataset")

# ── NanoString CosMx ─────────────────────────────────────────────────────────
print("\n--- NanoString CosMx ---")
pub("cosmx_nsclc_6k",
    "cosmx", "cosmx_6k", "spatial_transcriptomics", "probe_based",
    "https://nanostring.com/products/cosmx-spatial-molecular-imager/nsclc-data-release",
    "CosMx SMI human NSCLC; 6-cell-line dataset; 960 gene panel (public release)")

pub("cosmx_human_liver",
    "cosmx", "cosmx_6k", "spatial_transcriptomics", "probe_based",
    "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/liver-and-kidney",
    "CosMx SMI human liver + kidney; FFPE; 960 gene panel")

pub("cosmx_human_brain_ad",
    "cosmx", "cosmx_6k", "spatial_transcriptomics", "probe_based",
    "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/alzheimers-disease",
    "CosMx SMI human Alzheimers disease brain; 960 genes; FFPE")

pub("cosmx_human_gbm",
    "cosmx", "cosmx_6k", "spatial_transcriptomics", "probe_based",
    "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/glioblastoma",
    "CosMx SMI human glioblastoma; FFPE")

pub("cosmx_mouse_brain",
    "cosmx", "cosmx_6k", "spatial_transcriptomics", "probe_based",
    "https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/mouse-brain",
    "CosMx SMI mouse brain; FFPE; 960 gene panel")

# ── CODEX / MIBI / CyCIF ──────────────────────────────────────────────────────
print("\n--- CODEX / MIBI / CyCIF ---")
pub("codex_human_intestine",
    "codex", "codex", "spatial_proteomics", "multiplexed_fluorescence",
    "https://data.mendeley.com/datasets/mpjzbtfgfr/1",
    "Gut Cell Atlas CODEX; human intestine; 54 markers; Hickey et al. 2023 Cell")

pub("codex_human_tonsil",
    "codex", "codex", "spatial_proteomics", "multiplexed_fluorescence",
    "https://zenodo.org/record/5244551",
    "Goltsev et al. 2018 Cell; human tonsil CODEX; 27 markers")

pub("mibi_tnbc_tme",
    "mibi", "mibi", "spatial_proteomics", "multiplexed_fluorescence",
    "https://www.angelolab.com/mibi-data",
    "Angelo lab MIBI-TOF TNBC; 36 markers; Keren 2018 Cell")

pub("mibi_human_crc",
    "mibi", "mibi", "spatial_proteomics", "multiplexed_fluorescence",
    "https://zenodo.org/record/7626104",
    "MIBI human CRC; immune and tumor markers; Lin et al. 2023 Cancer Cell")

pub("cycif_human_crc",
    "cycif", "cycif", "imaging", "multiplexed_fluorescence",
    "https://www.synapse.org/#!Synapse:syn22345966",
    "Lin et al. 2023 Science; Human colorectal CyCIF atlas; 55 markers; HTAN CRC dataset")

pub("cycif_human_ovarian_cancer",
    "cycif", "cycif", "imaging", "multiplexed_fluorescence",
    "https://www.synapse.org/#!Synapse:syn17853080",
    "Labeling atlas of human ovarian cancer; CyCIF 60-plex")

pub("micsss_human_tumor",
    "micsss", "micsss", "imaging", "multiplexed_fluorescence",
    "https://zenodo.org/record/4890549",
    "MICSSS human colorectal tumor; multiple proteins; Eng et al. 2022")

# ── SlideSeq / Stereo-seq / Slide-tags ───────────────────────────────────────
print("\n--- SlideSeq / Stereo-seq / Slide-tags ---")
pub("slideseq_mouse_cerebellum",
    "slide_seq", "slide_seq", "spatial_transcriptomics", "sequencing_based",
    "https://singlecell.broadinstitute.org/single_cell/study/SCP354",
    "Slide-seq v1 mouse cerebellum; Rodriques et al. 2019 Science; 10um beads")

pub("slideseqv2_mouse_hippocampus",
    "slide_seq", "slide_seq", "spatial_transcriptomics", "sequencing_based",
    "https://singlecell.broadinstitute.org/single_cell/study/SCP815",
    "Slide-seq v2 mouse hippocampus; Stickels et al. 2021 Nat Biotechnol; ~10um")

pub("stereo_seq_mouse_embryo",
    "stereo_seq", "stereo_seq", "spatial_transcriptomics", "sequencing_based",
    "https://db.cngb.org/stomics/mosta",
    "Stereo-seq mouse embryo MOSTA; full developmental series; 500nm resolution")

pub("slide_tags_human_tonsil",
    "slide_tags", "slide_tags", "spatial_transcriptomics", "sequencing_based",
    "https://singlecell.broadinstitute.org/single_cell/study/SCP2169",
    "Slide-tags human tonsil; Russell et al. 2023 Nature; single-cell resolution")

# ── HTAN (Human Tumor Atlas Network) ─────────────────────────────────────────
print("\n--- HTAN ---")
pub("htan_crc_synapse",
    "various", "multimodal", "spatial_transcriptomics", "sequencing_based",
    "https://humantumoratlas.org",
    "HTAN human colorectal cancer atlas; Visium + CyCIF + scRNA; synapse.org/HTAN")

pub("htan_breast_cancer",
    "various", "multimodal", "spatial_transcriptomics", "sequencing_based",
    "https://humantumoratlas.org",
    "HTAN breast cancer atlas; multiple modalities; TCGA + spatial")

pub("htan_lung_cancer",
    "various", "multimodal", "spatial_transcriptomics", "sequencing_based",
    "https://humantumoratlas.org",
    "HTAN NSCLC; Visium + Xenium + scRNA")

pub("htan_melanoma",
    "various", "multimodal", "spatial_transcriptomics", "sequencing_based",
    "https://humantumoratlas.org",
    "HTAN melanoma; CODEX + scRNA; multi-center")

pub("htan_glioblastoma",
    "various", "multimodal", "spatial_transcriptomics", "sequencing_based",
    "https://humantumoratlas.org",
    "HTAN glioblastoma; spatial + scRNA")

pub("htan_ovarian_cancer",
    "various", "multimodal", "spatial_proteomics", "multiplexed_fluorescence",
    "https://humantumoratlas.org",
    "HTAN ovarian cancer; CyCIF + scRNA; OHSU center")

# ── HuBMAP ────────────────────────────────────────────────────────────────────
print("\n--- HuBMAP ---")
pub("hubmap_human_intestine",
    "various", "multimodal", "spatial_transcriptomics", "sequencing_based",
    "https://portal.hubmapconsortium.org",
    "HuBMAP human intestine; Visium + CODEX + scRNA; Gut Cell Atlas")

pub("hubmap_human_kidney",
    "various", "multimodal", "spatial_transcriptomics", "sequencing_based",
    "https://portal.hubmapconsortium.org",
    "HuBMAP human kidney; Visium + CODEX; Kidney Precision Medicine Project")

pub("hubmap_human_spleen",
    "various", "multimodal", "spatial_proteomics", "multiplexed_fluorescence",
    "https://portal.hubmapconsortium.org",
    "HuBMAP human spleen; CODEX + scRNA")

pub("hubmap_human_thymus",
    "various", "multimodal", "spatial_transcriptomics", "sequencing_based",
    "https://portal.hubmapconsortium.org",
    "HuBMAP human thymus developmental atlas; Visium + scRNA")

pub("hubmap_human_liver",
    "various", "multimodal", "spatial_transcriptomics", "sequencing_based",
    "https://portal.hubmapconsortium.org",
    "HuBMAP human liver; Visium + CODEX")

pub("hubmap_human_lung",
    "various", "multimodal", "spatial_transcriptomics", "sequencing_based",
    "https://portal.hubmapconsortium.org",
    "HuBMAP human lung; Visium + CODEX; LungMAP consortium")

pub("hubmap_human_heart",
    "various", "multimodal", "spatial_transcriptomics", "sequencing_based",
    "https://portal.hubmapconsortium.org",
    "HuBMAP human heart; Visium + CODEX; cardiac cell atlas")

# ── HEST-1k ────────────────────────────────────────────────────────────────────
print("\n--- HEST-1k ---")
pub("hest1k_collection",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    "https://huggingface.co/datasets/MahmoodLab/hest",
    "HEST-1k: 1108 spatial transcriptomics samples + H&E images; 131 datasets; huggingface MahmoodLab/hest")

pub("hest1k_cptac_breast",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    "https://huggingface.co/datasets/MahmoodLab/hest",
    "HEST-1k CPTAC breast cancer subset (>100 samples)")

pub("hest1k_tcga_brca",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    "https://huggingface.co/datasets/MahmoodLab/hest",
    "HEST-1k TCGA-BRCA Visium subset; matched H&E + spatial expression")

# ── SpatialData scverse ────────────────────────────────────────────────────────
print("\n--- SpatialData scverse tutorial datasets ---")
pub("spatialdata_mibitof_breast",
    "mibi", "mibi", "spatial_proteomics", "multiplexed_fluorescence",
    "https://s3.embl.de/spatialdata/spatialdata-sandbox/mibitof.zarr.zip",
    "SpatialData MIBI-TOF breast tutorial dataset (zarr format)")

pub("spatialdata_visium_hd_her2",
    "visium_hd", "visium_hd", "spatial_transcriptomics", "sequencing_based",
    "https://s3.embl.de/spatialdata/spatialdata-sandbox/visium_hd_3.0.0_io.zarr.zip",
    "SpatialData Visium HD HER2+ breast cancer tutorial (zarr format)")

pub("spatialdata_xenium_breast",
    "xenium", "xenium", "spatial_transcriptomics", "probe_based",
    "https://s3.embl.de/spatialdata/spatialdata-sandbox/xenium_rep1_io.zarr.zip",
    "SpatialData Xenium breast cancer rep1 tutorial (zarr format)")

pub("spatialdata_merfish_brain",
    "merfish", "merfish", "spatial_transcriptomics", "probe_based",
    "https://s3.embl.de/spatialdata/spatialdata-sandbox/merfish.zarr.zip",
    "SpatialData MERFISH mouse brain tutorial (zarr format)")

pub("spatialdata_cosmx_brain",
    "cosmx", "cosmx_6k", "spatial_transcriptomics", "probe_based",
    "https://s3.embl.de/spatialdata/spatialdata-sandbox/cosmx.zarr.zip",
    "SpatialData CosMx mouse brain tutorial (zarr format)")

pub("spatialdata_visium_hd_mouse_int",
    "visium_hd", "visium_hd", "spatial_transcriptomics", "sequencing_based",
    "https://s3.embl.de/spatialdata/spatialdata-sandbox/visium_hd_3.0.0_io_aligned.zarr.zip",
    "SpatialData Visium HD mouse intestine (zarr format)")

# ── GEO notable accessions ────────────────────────────────────────────────────
print("\n--- GEO notable accessions ---")
pub("geo_GSE144239_crc_visium",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144239",
    "GSE144239: Visium human CRC; Mayer et al. 2021; paired normal/tumor")

pub("geo_GSE176078_tnbc_visium",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078",
    "GSE176078: Wu et al. 2021 Nature Genetics; TNBC spatial + scRNA")

pub("geo_GSE203612_human_liver",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203612",
    "GSE203612: Visium human liver cancer; hepatocellular + intrahepatic CCA")

pub("geo_GSE166737_mouse_brain_spatiotemporal",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166737",
    "GSE166737: Mouse brain spatiotemporal transcriptomics; Visium")

pub("geo_GSE214989_xenium_breast_cancer",
    "xenium", "xenium", "spatial_transcriptomics", "probe_based",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214989",
    "GSE214989: First published Xenium breast cancer dataset; Janesick et al. 2023 Nat Commun")

pub("geo_GSE188729_imc_crc_atlas",
    "imc", "imc", "spatial_proteomics", "mass_spec_imaging",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188729",
    "GSE188729: Hartmann et al. IMC colorectal cancer atlas")

pub("geo_GSE149580_imc_covid19",
    "imc", "imc", "spatial_proteomics", "mass_spec_imaging",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149580",
    "GSE149580: Rendeiro et al. 2021 Nature; COVID-19 lung IMC")

pub("geo_GSE169249_human_heart_visium",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169249",
    "GSE169249: Kuppe et al. 2022 Nature; human cardiac scRNA + Visium")

pub("geo_GSE198353_merfish_hypothalamus",
    "merfish", "merfish", "spatial_transcriptomics", "probe_based",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198353",
    "GSE198353: MERFISH mouse hypothalamus; diet-perturbation")

pub("geo_GSE181771_visium_human_pcos",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181771",
    "GSE181771: Visium human endometrium; PCOS study")

# ── scProAtlas and other curated references ───────────────────────────────────
print("\n--- scProAtlas / reference atlases ---")
pub("scproatlas_human_CODEX",
    "codex", "codex", "spatial_proteomics", "multiplexed_fluorescence",
    "https://scproatlas.com",
    "scProAtlas: curated multiplexed protein spatial atlas (CODEX/MIBI/CyCIF)")

pub("human_cell_atlas_spatial",
    "various", "multimodal", "spatial_transcriptomics", "sequencing_based",
    "https://www.humancellatlas.org/biological-networks/spatial-biology",
    "Human Cell Atlas Spatial Biology collection; multiple organs and technologies")

pub("mouse_brain_atlas_zhang2023",
    "merfish", "merfish", "spatial_transcriptomics", "probe_based",
    "https://zenodo.org/record/8179024",
    "Zhang et al. 2023 Nature; mouse whole-brain 1100-gene MERFISH atlas")

pub("human_cortex_allen_10x",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    "https://portal.brain-map.org/atlases-and-data/rnaseq",
    "Allen Institute 10x Visium human cortex; multiple regions; scRNA + spatial")

pub("gut_cell_atlas",
    "various", "multimodal", "spatial_transcriptomics", "sequencing_based",
    "https://www.gutcellatlas.org",
    "Gut Cell Atlas; human intestine; scRNA + CODEX; Elmentaite et al. 2021")

pub("lung_cell_atlas",
    "various", "multimodal", "spatial_transcriptomics", "sequencing_based",
    "https://www.lungcellatlas.org",
    "Human Lung Cell Atlas; Travaglini et al. 2020 Nature; scRNA + Visium")

pub("heart_cell_atlas",
    "various", "multimodal", "spatial_transcriptomics", "sequencing_based",
    "https://www.heartcellatlas.org",
    "Human Heart Cell Atlas; Litvinukova et al. 2020 Nature; scRNA + Visium")

pub("kidney_cell_atlas",
    "various", "multimodal", "spatial_transcriptomics", "sequencing_based",
    "https://www.kidneycellatlas.org",
    "Human Kidney Cell Atlas; Stewart et al. 2019 Science; scRNA + spatial")

# ── Benchmark datasets ────────────────────────────────────────────────────────
print("\n--- Benchmark / methods evaluation datasets ---")
pub("benchmark_squidpy_imc_demo",
    "imc", "imc", "spatial_proteomics", "mass_spec_imaging",
    "https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_imc.html",
    "Squidpy IMC tutorial dataset; Schapiro lung/CRC; available via squidpy.datasets")

pub("benchmark_seurat_visium_brain",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    "https://satijalab.org/seurat/articles/spatial_vignette",
    "Seurat Visium mouse brain tutorial dataset; 10x mouse brain coronal section")

pub("benchmark_stereopy_merfish",
    "merfish", "merfish", "spatial_transcriptomics", "probe_based",
    "https://stereopy.readthedocs.io/en/latest/",
    "STOmics stereopy tutorial MERFISH dataset; mouse brain")

pub("benchmark_giotto_codex",
    "codex", "codex", "spatial_proteomics", "multiplexed_fluorescence",
    "https://giottosuite.readthedocs.io/en/master/subsections/datasets",
    "Giotto CODEX tutorial dataset; human CODEX intestine")

pub("benchmark_stlearn_breast",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    "https://stlearn.readthedocs.io/en/latest/tutorials",
    "stlearn Visium human breast cancer tutorial")

pub("benchmark_banksy_hippocampus",
    "slide_seq", "slide_seq", "spatial_transcriptomics", "sequencing_based",
    "https://github.com/prabhakarlab/Banksy_py",
    "BANKSY benchmarking data; mouse hippocampus Slide-seq v2; Singhal et al. 2024 Nat Genet")

pub("benchmark_paste_prostate",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    "https://github.com/raphael-group/paste",
    "PASTE spatial alignment benchmark; human prostate; 3 serial sections")

pub("benchmark_tangram_scRNA_spatial",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    "https://github.com/broadinstitute/Tangram",
    "Tangram deconvolution benchmark; Visium + scRNA paired datasets")

pub("benchmark_cell2location_lymph_node",
    "visium", "visium", "spatial_transcriptomics", "sequencing_based",
    "https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html",
    "Cell2location Visium human lymph node tutorial; Kleshchevnikov et al. 2022 Nat Biotechnol")

pub("benchmark_sparcle_rna_imc",
    "imc", "imc", "spatial_proteomics", "mass_spec_imaging",
    "https://github.com/dpeerlab/SPARCLE",
    "SPARCLE integration benchmark; paired RNA-seq + IMC; Goltsev et al.")

# ── Scanpy/scverse standard tutorial datasets ─────────────────────────────────
print("\n--- Scanpy / scverse standard datasets ---")
pub("scanpy_pbmc3k_tutorial",
    "10x_chromium", "scrna", "single_cell", "sequencing_based",
    "https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html",
    "Scanpy PBMC 3k tutorial dataset; Zheng et al. 2017 Nat Commun")

pub("scanpy_mouse_intestine",
    "10x_chromium", "scrna", "single_cell", "sequencing_based",
    "https://scanpy.readthedocs.io/en/stable/tutorials",
    "Scanpy mouse intestinal epithelium tutorial; 10x 3prime v3")

pub("scanpy_human_pancreas_integration",
    "10x_chromium", "scrna", "single_cell", "sequencing_based",
    "https://scanpy-tutorials.readthedocs.io/en/latest/integrating-data-using-ingest.html",
    "Scanpy human pancreas scRNA integration benchmark; multiple labs/protocols")

pub("scvi_cortex_tutorial",
    "10x_chromium", "scrna", "single_cell", "sequencing_based",
    "https://scvi-tools.org/en/stable/tutorials",
    "scVI mouse cortex tutorial dataset; Tasic et al. 2018 Nature")

pub("scvi_pbmc_covid",
    "10x_chromium", "scrna", "single_cell", "sequencing_based",
    "https://scvi-tools.org/en/stable/tutorials/notebooks/scrna/scarches_scvi_tools.html",
    "scVI scArches PBMC + COVID reference mapping tutorial")


# ──────────────────────────────────────────────────────────────────────────────
# Final summary
# ──────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("DONE. Registry summary:")
projs = r.list_projects()
internal = [p for p in projs if p.get("project_type") == "internal"]
external = [p for p in projs if p.get("project_type") == "external"]
print(f"  Total projects : {len(projs)}")
print(f"  Internal       : {len(internal)}")
print(f"  External/public: {len(external)}")
print("=" * 70)
