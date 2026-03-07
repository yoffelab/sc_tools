#!/usr/bin/env python
"""One-time migration: move external public projects → data_sources catalog.
Also registers all /athena/project-saha/ data and creates the IBD project.

Run once from repo root:
    python scripts/migrate_external_to_data_sources.py
"""

from __future__ import annotations

from sc_tools.registry import Registry

r = Registry()


def ds(
    name: str,
    uri: str,
    *,
    description: str | None = None,
    platform: str,
    domain: str,
    imaging_modality: str,
    source_type: str,
    organism: str = "human",
    tissue: str | None = None,
    disease: str | None = None,
    n_samples: int | None = None,
    publication: str | None = None,
    access_notes: str | None = None,
    status: str = "available",
) -> int:
    """Wrapper with print."""
    src_id = r.register_data_source(
        name, uri,
        description=description,
        platform=platform,
        domain=domain,
        imaging_modality=imaging_modality,
        source_type=source_type,
        organism=organism,
        tissue=tissue,
        disease=disease,
        n_samples=n_samples,
        publication=publication,
        access_notes=access_notes,
        status=status,
    )
    print(f"  [+] {name}")
    return src_id


# ==============================================================================
# 0. Delete the 142 external "projects" that should be data_sources
# ==============================================================================
print("\n=== Removing external projects (will be re-added as data_sources) ===")
with r._session() as sess:
    external = sess.query(r._Project).filter_by(project_type="external").all()
    names = [p.name for p in external]
    print(f"  Deleting {len(names)} external project rows…")
    for p in external:
        sess.delete(p)
    sess.commit()
print(f"  Deleted {len(names)} rows.")


# ==============================================================================
# 1. Collaborator data: /athena/project-saha/ on cayuga
# ==============================================================================
print("\n=== project-saha (cayuga) ===")

# The multi-tissue CosMx/Xenium study (34 samples in data/)
# Naming: SAHA_{modality}_{tissue}_{slide}_{batch}
# CP/CR = CosMx (RDS Seurat objects); XR = Xenium (RDS Seurat objects)

ds("saha_cosmx_multitissue",
   "cayuga:/athena/project-saha/data/",
   description="Multi-tissue CosMx 1k study; 34 samples; APE/COL/ILE/LN/PANC/STO/PRCA/PDAC; "
               "RDS Seurat objects per sample (flatFiles + RawFiles + TileDB + H&E TIFF)",
   platform="cosmx",
   domain="spatial_transcriptomics",
   imaging_modality="probe_based",
   source_type="hpc_collaborator",
   organism="human",
   tissue="colon,ileum,appendix,pancreas,lymph_node,stomach,prostate",
   disease="cancer,IBD",
   n_samples=34,
   access_notes="cayuga:/athena/project-saha/data/; owner jip2007@med.cornell.edu; "
                "Seurat RDS format — needs rpy2 conversion to AnnData",
)

ds("saha_cosmx_1k_16",
   "cayuga:/athena/project-saha/CosMx_1k_16/",
   description="CosMx 1k panel; 16 QC-normalized Seurat RDS objects",
   platform="cosmx",
   domain="spatial_transcriptomics",
   imaging_modality="probe_based",
   source_type="hpc_collaborator",
   organism="human",
   n_samples=16,
   access_notes="cayuga:/athena/project-saha/CosMx_1k_16/; Seurat RDS format",
)

ds("saha_cosmx_6k_4",
   "cayuga:/athena/project-saha/CosMx_6k_4/",
   description="CosMx 6k panel; 4 QC-normalized Seurat RDS objects",
   platform="cosmx",
   domain="spatial_transcriptomics",
   imaging_modality="probe_based",
   source_type="hpc_collaborator",
   organism="human",
   n_samples=4,
   access_notes="cayuga:/athena/project-saha/CosMx_6k_4/; Seurat RDS format",
)

ds("saha_xenium_5k_4",
   "cayuga:/athena/project-saha/Xenium_5K_4/",
   description="Xenium 5k-gene panel; 4 QC-normalized Seurat RDS objects",
   platform="xenium",
   domain="spatial_transcriptomics",
   imaging_modality="probe_based",
   source_type="hpc_collaborator",
   organism="human",
   n_samples=4,
   access_notes="cayuga:/athena/project-saha/Xenium_5K_4/",
)

ds("saha_xenium_colon_4",
   "cayuga:/athena/project-saha/Xenium_colon_4/",
   description="Xenium colon; 4 samples with cell segmentation; Seurat RDS",
   platform="xenium",
   domain="spatial_transcriptomics",
   imaging_modality="probe_based",
   source_type="hpc_collaborator",
   organism="human",
   tissue="colon",
   n_samples=4,
   access_notes="cayuga:/athena/project-saha/Xenium_colon_4/",
)

ds("saha_xenium_mt_16",
   "cayuga:/athena/project-saha/Xenium_MT_16/",
   description="Xenium multi-tissue; 16 samples; Seurat RDS",
   platform="xenium",
   domain="spatial_transcriptomics",
   imaging_modality="probe_based",
   source_type="hpc_collaborator",
   organism="human",
   tissue="colon,ileum,appendix,pancreas",
   n_samples=16,
   access_notes="cayuga:/athena/project-saha/Xenium_MT_16/",
)

ds("saha_xenium_mt_noseg_4",
   "cayuga:/athena/project-saha/Xenium_MT_noseg_4/",
   description="Xenium MT 377-gene panel; 4 samples; no cell segmentation; Seurat RDS",
   platform="xenium",
   domain="spatial_transcriptomics",
   imaging_modality="probe_based",
   source_type="hpc_collaborator",
   organism="human",
   n_samples=4,
   access_notes="cayuga:/athena/project-saha/Xenium_MT_noseg_4/",
)

ds("saha_xenium_mt_withseg_4",
   "cayuga:/athena/project-saha/Xenium_MT_withseg_4/",
   description="Xenium MT 377-gene panel; 4 samples; with cell segmentation; Seurat RDS",
   platform="xenium",
   domain="spatial_transcriptomics",
   imaging_modality="probe_based",
   source_type="hpc_collaborator",
   organism="human",
   n_samples=4,
   access_notes="cayuga:/athena/project-saha/Xenium_MT_withseg_4/",
)

ds("saha_he_tif",
   "cayuga:/athena/project-saha/HE/",
   description="H&E TIFF images; companion to CosMx/Xenium samples (I0262, I0275-I0294…)",
   platform="he",
   domain="imaging",
   imaging_modality="brightfield",
   source_type="hpc_collaborator",
   organism="human",
   access_notes="cayuga:/athena/project-saha/HE/; TIFF format",
)

ds("saha_he_fullres_ndpi",
   "cayuga:/athena/project-saha/HE_fullres/",
   description="Full-resolution H&E NDPI scans (Hamamatsu format)",
   platform="he",
   domain="imaging",
   imaging_modality="brightfield",
   source_type="hpc_collaborator",
   organism="human",
   access_notes="cayuga:/athena/project-saha/HE_fullres/; NDPI format — needs QuPath/openslide",
)

ds("saha_ibd_scrna",
   "cayuga:/athena/project-saha/SAHA_IBD_RNA.h5ad",
   description="IBD scRNA-seq reference; 11 GB h5ad; multi-sample; "
               "likely used as deconvolution reference for IBD spatial project",
   platform="scrna",
   domain="single_cell",
   imaging_modality="sequencing_based",
   source_type="hpc_collaborator",
   organism="human",
   tissue="colon,ileum",
   disease="IBD",
   access_notes="cayuga:/athena/project-saha/SAHA_IBD_RNA.h5ad; h5ad format; 11 GB",
)

ds("saha_ibd_spatial_pending",
   "cayuga:/athena/project-saha/data_IBD/",
   description="IBD spatial data directory (currently empty; data expected from Saha lab)",
   platform="cosmx",
   domain="spatial_transcriptomics",
   imaging_modality="probe_based",
   source_type="hpc_collaborator",
   organism="human",
   tissue="colon,ileum",
   disease="IBD",
   access_notes="cayuga:/athena/project-saha/data_IBD/; directory created 2026-03-06; "
                "contact jip2007 for data delivery timeline",
   status="pending_download",
)

ds("saha_metadata_nc_compare",
   "cayuga:/athena/project-saha/NC_compare_11122025_meta.csv",
   description="Sample metadata / NC comparison table (dated 2025-11-12)",
   platform="metadata",
   domain="single_cell",
   imaging_modality="sequencing_based",
   source_type="hpc_collaborator",
   organism="human",
   tissue="colon,ileum,appendix,pancreas",
   disease="IBD,cancer",
   access_notes="cayuga:/athena/project-saha/NC_compare_11122025_meta.csv",
)


# ==============================================================================
# 2. Create the IBD spatial project and link its data sources
# ==============================================================================
print("\n=== Creating IBD spatial project ===")
r.add_project(
    "ibd_spatial",
    platform="cosmx",
    data_type="cosmx",
    domain="spatial_transcriptomics",
    imaging_modality="probe_based",
    project_type="external",   # collaboration with Saha lab
    visibility="private",
)
print("  [+] project: ibd_spatial")

# Link data sources to the project
for src_name, role in [
    ("saha_ibd_spatial_pending",  "input"),
    ("saha_cosmx_multitissue",    "input"),
    ("saha_cosmx_1k_16",          "input"),
    ("saha_cosmx_6k_4",           "input"),
    ("saha_xenium_colon_4",       "input"),
    ("saha_xenium_mt_16",         "input"),
    ("saha_ibd_scrna",            "reference"),
    ("saha_metadata_nc_compare",  "supplementary"),
    ("saha_he_tif",               "supplementary"),
]:
    r.link_project_data_source("ibd_spatial", src_name, role=role)
    print(f"  [~] linked: ibd_spatial → {src_name} ({role})")


# ==============================================================================
# 3. Register key public data sources (no longer as projects)
# ==============================================================================
print("\n=== Public data sources ===")

# ── 10x Genomics ─────────────────────────────────────────────────────────────
BASE = "https://www.10xgenomics.com/datasets"

for name, uri_tail, pl, tissue, disease, n in [
    ("10x_visium_human_breast_cancer",      "/spatial-gene-expression-v1-human-breast-cancer-1",          "visium",    "breast",    "breast_cancer", 2),
    ("10x_visium_human_heart",              "/spatial-gene-expression-human-heart",                        "visium",    "heart",     None,            1),
    ("10x_visium_human_brain_coronal",      "/spatial-gene-expression-human-cerebral-cortex",              "visium",    "brain",     None,            1),
    ("10x_visium_human_lymph_node",         "/spatial-gene-expression-human-lymph-node",                   "visium",    "lymph_node",None,            1),
    ("10x_visium_human_ovarian_cancer",     "/spatial-gene-expression-v1-human-ovarian-cancer",            "visium",    "ovary",     "ovarian_cancer",1),
    ("10x_visium_mouse_brain_coronal",      "/spatial-gene-expression-mouse-brain-coronal-section-1",      "visium",    "brain",     None,            1),
    ("10x_visium_human_colon_cancer",       "/spatial-gene-expression-human-colon-cancer-1",               "visium",    "colon",     "colon_cancer",  1),
    ("10x_visium_human_glioblastoma",       "/spatial-gene-expression-human-glioblastoma",                 "visium",    "brain",     "glioblastoma",  1),
    ("10x_visium_human_prostate_cancer",    "/spatial-gene-expression-human-prostate-cancer",              "visium",    "prostate",  "prostate_cancer",1),
    ("10x_visium_human_skin",               "/spatial-gene-expression-human-skin",                         "visium",    "skin",      None,            1),
    ("10x_visium_cytassist_human_lung",     "/visium-cytassist-gene-expression-human-lung-cancer",         "visium",    "lung",      "lung_cancer",   1),
    ("10x_visium_cytassist_human_crc",      "/visium-cytassist-gene-expression-human-colorectal-cancer",   "visium",    "colon",     "colorectal_cancer",1),
    ("10x_visium_cytassist_mouse_brain",    "/visium-cytassist-gene-expression-mouse-brain",               "visium",    "brain",     None,            1),
    ("10x_visium_cytassist_human_melanoma", "/visium-cytassist-gene-expression-human-melanoma",            "visium",    "skin",      "melanoma",      1),
    ("10x_visium_cytassist_human_tonsil",   "/visium-cytassist-gene-expression-human-tonsil",              "visium",    "tonsil",    None,            1),
    ("10x_visium_cytassist_human_kidney",   "/visium-cytassist-gene-expression-human-kidney",              "visium",    "kidney",    None,            1),
    ("10x_visium_ffpe_human_breast",        "/spatial-gene-expression-v1-human-breast-cancer-ffpe",        "visium",    "breast",    "breast_cancer", 1),
    ("10x_visium_ffpe_human_cerebellum",    "/spatial-gene-expression-human-cerebellum-ffpe",              "visium",    "brain",     None,            1),
    ("10x_visium_hd_human_crc",             "/visium-hd-cytassist-gene-expression-libraries-of-human-crc","visium_hd", "colon",     "colorectal_cancer",1),
    ("10x_visium_hd_human_breast_cancer",   "/visium-hd-cytassist-gene-expression-libraries-of-human-breast-cancer","visium_hd","breast","breast_cancer",1),
    ("10x_visium_hd_mouse_brain",           "/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain","visium_hd","brain",    None,            1),
    ("10x_visium_hd_mouse_intestine",       "/visium-hd-cytassist-gene-expression-libraries-of-mouse-intestine","visium_hd","intestine",None,         1),
    ("10x_visium_hd_human_pancreatic_cancer","/visium-hd-cytassist-gene-expression-libraries-of-human-pancreatic-cancer","visium_hd","pancreas","pancreatic_cancer",1),
    ("10x_visium_hd_human_lung_cancer",     "/visium-hd-cytassist-gene-expression-libraries-of-human-lung-cancer","visium_hd","lung","lung_cancer",  1),
    ("10x_visium_hd_post_xenium_human_bc",  "/visium-hd-cytassist-gene-expression-libraries-of-post-xenium-human-breast-cancer","visium_hd","breast","breast_cancer",1),
    ("10x_xenium_human_breast_cancer",      "/xenium-human-breast-cancer-gene-expression",                 "xenium",    "breast",    "breast_cancer", 2),
    ("10x_xenium_human_lung_cancer",        "/xenium-human-lung-cancer-gene-expression",                   "xenium",    "lung",      "lung_cancer",   1),
    ("10x_xenium_human_brain",              "/xenium-human-brain-gene-expression",                         "xenium",    "brain",     None,            1),
    ("10x_xenium_human_colon_cancer",       "/xenium-human-colon-cancer-gene-expression",                  "xenium",    "colon",     "colon_cancer",  1),
    ("10x_xenium_human_pancreatic_cancer",  "/xenium-human-pancreatic-cancer-gene-expression",             "xenium",    "pancreas",  "pancreatic_cancer",1),
    ("10x_xenium_human_lymph_node",         "/xenium-human-lymph-node-gene-expression",                    "xenium",    "lymph_node",None,            1),
    ("10x_xenium_mouse_brain",              "/xenium-mouse-brain-gene-expression",                         "xenium",    "brain",     None,            1),
    ("10x_xenium_human_crc_5k",             "/xenium-human-colorectal-cancer-gene-expression-5000-panel",  "xenium",    "colon",     "colorectal_cancer",1),
    ("10x_xenium_human_melanoma",           "/xenium-human-melanoma-gene-expression",                      "xenium",    "skin",      "melanoma",      1),
    ("10x_pbmc_3k",                         "/pbmc-3k-filtered-feature-bc-matrix",                         "10x_chromium","blood",  None,            1),
    ("10x_pbmc_10k_v3",                     "/10k-human-pbmcs-3-v3-1-chromium-x",                          "10x_chromium","blood",  None,            1),
]:
    org = "mouse" if "mouse" in name else "human"
    ds(name, BASE + uri_tail,
       platform=pl,
       domain="single_cell" if pl == "10x_chromium" else "spatial_transcriptomics",
       imaging_modality="sequencing_based" if pl in ("visium","visium_hd","10x_chromium") else "probe_based",
       source_type="public_10x",
       organism=org,
       tissue=tissue,
       disease=disease,
       n_samples=n,
       status="pending_download",
    )

# ── IMC public (Zenodo / imcdatasets) ─────────────────────────────────────────
for name, uri, tissue, disease, pub in [
    ("imc_damond2019_pancreas",   "https://zenodo.org/record/3951011",  "pancreas","T1D","Damond et al. 2019 Cell Metab doi:10.1016/j.cmet.2019.01.021"),
    ("imc_jackson2020_breast",    "https://zenodo.org/record/7575859",  "breast",  "breast_cancer","Jackson & Fischer 2020 Nature doi:10.1038/s41586-019-1876-x"),
    ("imc_schapiro2022_lung_crc", "https://zenodo.org/record/5949116",  "lung,colon","cancer","Schapiro et al. 2022 Nat Methods doi:10.1038/s41592-021-01318-y"),
    ("imc_hoch2022_melanoma",     "https://zenodo.org/record/7637156",  "skin",    "melanoma","Hoch et al. 2022 Sci Immunol doi:10.1126/sciimmunol.abk0892"),
    ("imc_crc_atlas",             "https://zenodo.org/record/6024273",  "colon",   "colorectal_cancer","Hartmann et al. 2021 Cancer Cell doi:10.1016/j.ccell.2020.07.005"),
    ("imc_basel_zurich_breast",   "https://zenodo.org/record/3518284",  "breast",  "breast_cancer","Ali et al. 2020 Nat Cancer doi:10.1038/s43018-020-0026-6"),
    ("imc_covid19_lung",          "https://zenodo.org/record/5018260",  "lung",    "COVID-19","Rendeiro et al. 2021 Nature doi:10.1038/s41586-021-03475-6"),
    ("imc_keren2018_tnbc",        "https://www.angelolab.com/mibi-data","breast",  "TNBC","Keren et al. 2018 Cell doi:10.1016/j.cell.2018.08.039"),
]:
    ds(name, uri,
       platform="imc", domain="spatial_proteomics", imaging_modality="mass_spec_imaging",
       source_type="public_zenodo",
       organism="human", tissue=tissue, disease=disease,
       publication=pub, status="pending_download",
    )

# ── MERFISH / Vizgen ──────────────────────────────────────────────────────────
for name, uri, platform, tissue, disease, st in [
    ("merfish_mouse_brain_abc_atlas", "https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas",
     "merfish","brain",None,"public_portal"),
    ("merfish_mouse_hypothalamus",    "https://datadryad.org/stash/dataset/doi:10.5061/dryad.8t8s248",
     "merfish","brain",None,"public_portal"),
    ("vizgen_merscope_human_liver",   "https://console.cloud.google.com/storage/browser/vz-liver-showcase",
     "merscope","liver","liver_cancer","public_portal"),
    ("vizgen_merscope_mouse_brain",   "https://console.cloud.google.com/storage/browser/vz-brain-showcase",
     "merscope","brain",None,"public_portal"),
]:
    org = "mouse" if "mouse" in name else "human"
    ds(name, uri,
       platform=platform, domain="spatial_transcriptomics", imaging_modality="probe_based",
       source_type=st, organism=org, tissue=tissue, disease=disease, status="pending_download",
    )

# ── CosMx public ─────────────────────────────────────────────────────────────
for name, uri, tissue, disease in [
    ("cosmx_nsclc_6k",   "https://nanostring.com/products/cosmx-spatial-molecular-imager/nsclc-data-release","lung","lung_cancer"),
    ("cosmx_human_brain_ad","https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/alzheimers-disease","brain","Alzheimers"),
]:
    ds(name, uri,
       platform="cosmx", domain="spatial_transcriptomics", imaging_modality="probe_based",
       source_type="public_portal", organism="human", tissue=tissue, disease=disease,
       status="pending_download",
    )

# ── CODEX / CyCIF ─────────────────────────────────────────────────────────────
for name, uri, pl, tissue, disease, pub_str in [
    ("codex_human_intestine","https://data.mendeley.com/datasets/mpjzbtfgfr/1","codex","intestine",None,"Hickey et al. 2023 Cell doi:10.1016/j.cell.2023.01.009"),
    ("mibi_tnbc_tme",        "https://www.angelolab.com/mibi-data","mibi","breast","TNBC","Keren et al. 2018 Cell doi:10.1016/j.cell.2018.08.039"),
    ("cycif_human_crc",      "https://www.synapse.org/#!Synapse:syn22345966","cycif","colon","colorectal_cancer","Lin et al. 2023 Science doi:10.1126/science.abi7357"),
]:
    ds(name, uri,
       platform=pl, domain="spatial_proteomics", imaging_modality="multiplexed_fluorescence",
       source_type="public_portal", organism="human", tissue=tissue, disease=disease,
       publication=pub_str, status="pending_download",
    )

# ── Portal collections ────────────────────────────────────────────────────────
for name, uri, dom, tissue, disease in [
    ("htan_crc",      "https://humantumoratlas.org","spatial_transcriptomics","colon","colorectal_cancer"),
    ("htan_breast",   "https://humantumoratlas.org","spatial_transcriptomics","breast","breast_cancer"),
    ("htan_lung",     "https://humantumoratlas.org","spatial_transcriptomics","lung","lung_cancer"),
    ("htan_melanoma", "https://humantumoratlas.org","imaging","skin","melanoma"),
    ("hubmap_intestine","https://portal.hubmapconsortium.org","spatial_transcriptomics","intestine",None),
    ("hubmap_kidney",   "https://portal.hubmapconsortium.org","spatial_transcriptomics","kidney",None),
    ("hubmap_lung",     "https://portal.hubmapconsortium.org","spatial_transcriptomics","lung",None),
    ("hest1k_collection","https://huggingface.co/datasets/MahmoodLab/hest","spatial_transcriptomics","various","various"),
    ("gut_cell_atlas",  "https://www.gutcellatlas.org","spatial_transcriptomics","intestine",None),
    ("lung_cell_atlas", "https://www.lungcellatlas.org","spatial_transcriptomics","lung",None),
]:
    ds(name, uri,
       platform="various", domain=dom, imaging_modality="sequencing_based",
       source_type="public_portal", organism="human", tissue=tissue, disease=disease,
       status="pending_download",
    )

# ── GEO accessions ────────────────────────────────────────────────────────────
for name, acc, pl, tissue, disease, pub_str in [
    ("geo_GSE144239_crc_visium",    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144239","visium","colon","colorectal_cancer","Mayer et al. 2021 GSE144239"),
    ("geo_GSE176078_tnbc_visium",   "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078","visium","breast","TNBC","Wu et al. 2021 Nat Genet doi:10.1038/s41588-021-00911-1"),
    ("geo_GSE214989_xenium_breast", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214989","xenium","breast","breast_cancer","Janesick et al. 2023 Nat Commun doi:10.1038/s41467-023-43458-5"),
    ("geo_GSE188729_imc_crc",       "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188729","imc","colon","colorectal_cancer","Hartmann et al. GSE188729"),
    ("geo_GSE169249_heart_visium",  "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169249","visium","heart","cardiac_fibrosis","Kuppe et al. 2022 Nature doi:10.1038/s41586-022-05060-x"),
    ("geo_GSE203612_liver_visium",  "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203612","visium","liver","liver_cancer","GSE203612"),
]:
    org = "mouse" if "mouse" in name else "human"
    ds(name, acc,
       platform=pl,
       domain="spatial_proteomics" if pl == "imc" else "spatial_transcriptomics",
       imaging_modality="mass_spec_imaging" if pl == "imc" else ("probe_based" if pl == "xenium" else "sequencing_based"),
       source_type="public_geo",
       organism=org, tissue=tissue, disease=disease, publication=pub_str,
       status="pending_download",
    )

# ── SpatialData / scverse sandbox ────────────────────────────────────────────
for name, uri, pl, tissue in [
    ("spatialdata_xenium_breast",   "https://s3.embl.de/spatialdata/spatialdata-sandbox/xenium_rep1_io.zarr.zip","xenium","breast"),
    ("spatialdata_visium_hd_her2",  "https://s3.embl.de/spatialdata/spatialdata-sandbox/visium_hd_3.0.0_io.zarr.zip","visium_hd","breast"),
    ("spatialdata_merfish_brain",   "https://s3.embl.de/spatialdata/spatialdata-sandbox/merfish.zarr.zip","merfish","brain"),
    ("spatialdata_mibitof_breast",  "https://s3.embl.de/spatialdata/spatialdata-sandbox/mibitof.zarr.zip","mibi","breast"),
]:
    org = "mouse" if "mouse" in name else "human"
    ds(name, uri,
       platform=pl, domain="spatial_transcriptomics", imaging_modality="probe_based" if pl in ("xenium","merfish") else "sequencing_based",
       source_type="public_portal", organism=org, tissue=tissue,
       status="pending_download",
    )

# ── scRNA reference datasets ──────────────────────────────────────────────────
for name, uri, tissue in [
    ("scanpy_pbmc3k",    "https://scanpy.readthedocs.io/en/stable/tutorials","blood"),
    ("scvi_cortex",      "https://scvi-tools.org/en/stable/tutorials","brain"),
]:
    ds(name, uri,
       platform="10x_chromium", domain="single_cell", imaging_modality="sequencing_based",
       source_type="public_portal", organism="human", tissue=tissue,
       status="pending_download",
    )

# ==============================================================================
# Final summary
# ==============================================================================
print("\n=== Final counts ===")
from sc_tools.registry import Registry as _R
r2 = _R()
projs = r2.list_projects()
sources = r2.list_data_sources()
int_p = [p for p in projs if p.get("project_type") == "internal"]
ext_p = [p for p in projs if p.get("project_type") == "external"]
hpc_s = [s for s in sources if "hpc" in (s.get("source_type") or "")]
pub_s = [s for s in sources if "public" in (s.get("source_type") or "")]

print(f"  Projects  — internal: {len(int_p)}, external/collab: {len(ext_p)}")
print(f"  DataSources — HPC: {len(hpc_s)}, public: {len(pub_s)}, total: {len(sources)}")
print("\nDone.")
