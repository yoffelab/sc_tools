# BioData Type Hierarchy

The platform registry in `sc_tools.biodata` organizes biological data platforms into a three-level hierarchy:

```
BioDataType (5 JTI types)
  -> Modality (human-readable family)
    -> Platform (slug + version)
```

- **BioDataType** maps to a SQLAlchemy Joined Table Inheritance (JTI) subclass in the registry ORM (`bio_data` base table + type-specific child table).
- **Modality** is a human-readable grouping string (e.g. `Spatial Proteomics - Mass Spec`) shared by platforms with similar experimental design.
- **Platform** is a unique slug (e.g. `imc`, `xenium`, `chromium_3p`) with vendor, resolution, measurement type, and default field values stored in a `PlatformSpec` dataclass.

## ORM Subtype Mapping

| `biodata_type` | ORM class | Child table | Measurement |
|---|---|---|---|
| `spatial_seq` | `SpatialSeqData` | `spatial_seq_data` | RNA (spatial) |
| `image` | `BioImage` | `bio_images` | Protein / morphology |
| `rnaseq` | `RNASeqData` | `rnaseq_data` | RNA (dissociated / bulk) |
| `epigenomics` | `EpigenomicsData` | `epigenomics_data` | Chromatin |
| `genome_seq` | `GenomeSeqData` | `genome_seq_data` | DNA |

---

## Complete Modality Taxonomy

| `biodata_type` | Modality | Example platforms |
|---|---|---|
| `spatial_seq` | Spatial Transcriptomics - Sequencing | `visium`, `visium_hd`, `stereo_seq`, `slide_seq` |
| `spatial_seq` | Spatial Transcriptomics - ISH | `xenium`, `cosmx_1k`, `merscope`, `seqfish` |
| `spatial_seq` | Spatial Transcriptomics - Region Capture | `geomx`, `visium_cytassist` |
| `image` | Spatial Proteomics - Mass Spec | `imc`, `mibi`, `maldi_ims` |
| `image` | Spatial Proteomics - Cyclic IF | `phenocycler`, `macsima`, `cycif`, `comet`, `orion` |
| `image` | Spatial Proteomics - Standard IF/IHC | `mihc`, `if_standard`, `vectra_polaris` |
| `image` | Histology | `he`, `brightfield`, `phase_contrast`, `pas`, `masson_trichrome`, `electron_microscopy` |
| `rnaseq` | Single-Cell RNA-seq - Droplet | `chromium_3p`, `chromium_5p`, `indrops`, `dropseq` |
| `rnaseq` | Single-Cell RNA-seq - Combinatorial | `parse_evercode`, `scale_rna`, `split_seq` |
| `rnaseq` | Single-Cell RNA-seq - Microwell | `bd_rhapsody`, `seq_well`, `microwell_seq` |
| `rnaseq` | Single-Cell RNA-seq - Plate | `smart_seq2`, `smart_seq3`, `cel_seq2` |
| `rnaseq` | Bulk RNA-seq | `illumina_rnaseq`, `bgi_rnaseq`, `ont_direct_rna`, `pacbio_isoseq` |
| `rnaseq` | Multiomics | `cite_seq`, `dogma_seq`, `share_seq`, `asap_seq`, `tea_seq` |
| `epigenomics` | Epigenomics - Bulk | `atac_seq`, `chip_seq`, `cut_and_run`, `bisulfite_seq`, `hi_c` |
| `epigenomics` | Epigenomics - Single Cell | `sc_atac_seq`, `sc_cut_and_tag`, `sc_methylation` |
| `epigenomics` | Epigenomics - Spatial | `spatial_atac`, `spatial_cut_tag` |
| `genome_seq` | Genome Sequencing - Short Read | `illumina_wgs`, `illumina_wes`, `bgi_wgs` |
| `genome_seq` | Genome Sequencing - Long Read | `pacbio_hifi`, `ont_wgs`, `element_aviti` |
| `genome_seq` | Genome Sequencing - Targeted | `targeted_panel`, `amplicon_seq` |

---

## Platform Listing by Modality

### Spatial Transcriptomics - Sequencing (`spatial_seq`)

| Slug | Label | Resolution | Vendor |
|---|---|---|---|
| `visium` | 10x Visium | spot | 10x_genomics |
| `visium_hd` | 10x Visium HD | spot | 10x_genomics |
| `visium_hd_cell` | 10x Visium HD Cell | single_cell | 10x_genomics |
| `stereo_seq` | BGI Stereo-seq | subcellular | bgi |
| `slide_seq` | Slide-seq / Slide-seqV2 | spot | broad_institute |
| `hdst` | High-Definition Spatial Transcriptomics | spot | broad_institute |
| `pixel_seq` | PIXEL-seq | subcellular | academic |
| `seq_scope` | Seq-Scope | subcellular | academic |
| `open_st` | Open-ST | subcellular | academic |
| `dbit_seq` | DBiT-seq | spot | academic |
| `sci_space` | sci-Space | spot | academic |
| `magic_seq` | MAGIC-seq | spot | academic |

### Spatial Transcriptomics - ISH (`spatial_seq`)

| Slug | Label | Resolution | Vendor |
|---|---|---|---|
| `xenium` | 10x Xenium | subcellular | 10x_genomics |
| `cosmx_1k` | NanoString CosMx SMI 1k | single_cell | nanostring |
| `cosmx_6k` | NanoString CosMx SMI 6k | single_cell | nanostring |
| `cosmx_wt` | NanoString CosMx Whole Transcriptome | single_cell | nanostring |
| `merscope` | Vizgen MERSCOPE / MERFISH | subcellular | vizgen |
| `seqfish` | seqFISH / seqFISH+ | subcellular | academic |
| `starmap` | STARmap / STARmap PLUS | subcellular | academic |
| `hybiss` | HybISS | subcellular | academic |
| `eelfish` | EEL FISH | subcellular | academic |
| `resolve` | Resolve Biosciences Molecular Cartography | subcellular | resolve |

### Spatial Transcriptomics - Region Capture (`spatial_seq`)

| Slug | Label | Resolution | Vendor |
|---|---|---|---|
| `geomx` | NanoString GeoMx DSP | spot | nanostring |
| `visium_cytassist` | 10x CytAssist | spot | 10x_genomics |

### Spatial Proteomics - Mass Spec (`image`)

| Slug | Label | Resolution | Vendor |
|---|---|---|---|
| `imc` | Imaging Mass Cytometry | single_cell | standard_biotools |
| `mibi` | Multiplexed Ion Beam Imaging | single_cell | ionpath |
| `maldi_ims` | MALDI Imaging Mass Spec | single_cell | bruker |

### Spatial Proteomics - Cyclic IF (`image`)

| Slug | Label | Resolution | Vendor |
|---|---|---|---|
| `phenocycler` | Akoya PhenoCycler (CODEX) | single_cell | akoya |
| `macsima` | Miltenyi MACSima | single_cell | miltenyi |
| `cycif` | CyCIF | single_cell | academic |
| `insituplex` | Ultivue InSituPlex | single_cell | ultivue |
| `comet` | Lunaphore COMET | single_cell | lunaphore |
| `orion` | RareCyte ORION | single_cell | rarecyte |

### Spatial Proteomics - Standard IF/IHC (`image`)

| Slug | Label | Resolution | Vendor |
|---|---|---|---|
| `mihc` | Multiplex IHC (Vectra/Polaris) | single_cell | akoya |
| `if_standard` | Standard Immunofluorescence | single_cell | generic |
| `vectra_polaris` | Akoya Vectra Polaris | single_cell | akoya |

### Histology (`image`)

| Slug | Label | Resolution | Vendor |
|---|---|---|---|
| `he` | H&E Staining | bulk | generic |
| `brightfield` | Unstained Brightfield | bulk | generic |
| `phase_contrast` | Phase Contrast Microscopy | bulk | generic |
| `pas` | PAS Staining | bulk | generic |
| `masson_trichrome` | Masson Trichrome | bulk | generic |
| `electron_microscopy` | Electron Microscopy (TEM/SEM) | bulk | generic |

### Single-Cell RNA-seq - Droplet (`rnaseq`)

| Slug | Label | Resolution | Vendor |
|---|---|---|---|
| `chromium_3p` | 10x Chromium 3' | single_cell | 10x_genomics |
| `chromium_5p` | 10x Chromium 5' | single_cell | 10x_genomics |
| `chromium_flex` | 10x Chromium Flex | single_cell | 10x_genomics |
| `chromium_multiome` | 10x Multiome (RNA + ATAC) | single_cell | 10x_genomics |
| `indrops` | inDrops | single_cell | academic |
| `dropseq` | Drop-seq | single_cell | academic |

### Single-Cell RNA-seq - Combinatorial (`rnaseq`)

| Slug | Label | Resolution | Vendor |
|---|---|---|---|
| `parse_evercode` | Parse Biosciences Evercode | single_cell | parse |
| `scale_rna` | Scale Biosciences | single_cell | scale |
| `split_seq` | SPLiT-seq | single_cell | academic |

### Single-Cell RNA-seq - Microwell (`rnaseq`)

| Slug | Label | Resolution | Vendor |
|---|---|---|---|
| `bd_rhapsody` | BD Rhapsody | single_cell | bd |
| `seq_well` | Seq-Well | single_cell | academic |
| `microwell_seq` | Microwell-seq | single_cell | academic |

### Single-Cell RNA-seq - Plate (`rnaseq`)

| Slug | Label | Resolution | Vendor |
|---|---|---|---|
| `smart_seq2` | Smart-seq2 | single_cell | academic |
| `smart_seq3` | Smart-seq3 | single_cell | academic |
| `cel_seq2` | CEL-Seq2 | single_cell | academic |

### Bulk RNA-seq (`rnaseq`)

| Slug | Label | Resolution | Vendor |
|---|---|---|---|
| `illumina_rnaseq` | Illumina Short-Read RNA-seq | bulk | illumina |
| `bgi_rnaseq` | BGI DNBSEQ RNA-seq | bulk | bgi |
| `ont_direct_rna` | Oxford Nanopore Direct RNA-seq | bulk | ont |
| `pacbio_isoseq` | PacBio Iso-Seq | bulk | pacbio |
| `ultima_rnaseq` | Ultima Genomics UG 100 | bulk | ultima |

### Multiomics (`rnaseq`)

| Slug | Label | Resolution | Vendor |
|---|---|---|---|
| `cite_seq` | CITE-seq | single_cell | biolegend |
| `reap_seq` | REAP-seq | single_cell | academic |
| `asap_seq` | ASAP-seq | single_cell | academic |
| `share_seq` | SHARE-seq | single_cell | academic |
| `dogma_seq` | DOGMA-seq | single_cell | academic |
| `tea_seq` | TEA-seq | single_cell | academic |
| `snare_seq` | SNARE-seq | single_cell | academic |
| `paired_seq` | Paired-seq | single_cell | academic |

### Epigenomics - Bulk (`epigenomics`)

| Slug | Label | Resolution | Vendor |
|---|---|---|---|
| `atac_seq` | ATAC-seq | bulk | generic |
| `chip_seq` | ChIP-seq | bulk | generic |
| `cut_and_run` | CUT&RUN | bulk | generic |
| `cut_and_tag` | CUT&Tag | bulk | generic |
| `bisulfite_seq` | Bisulfite Sequencing (WGBS/RRBS) | bulk | generic |
| `em_seq` | EM-seq | bulk | generic |
| `hi_c` | Hi-C | bulk | generic |

### Epigenomics - Single Cell (`epigenomics`)

| Slug | Label | Resolution | Vendor |
|---|---|---|---|
| `sc_atac_seq` | scATAC-seq | single_cell | generic |
| `sc_cut_and_tag` | Single-cell CUT&Tag | single_cell | generic |
| `sc_methylation` | Single-cell Bisulfite (scBS-seq) | single_cell | generic |
| `sc_hi_c` | Single-cell Hi-C | single_cell | generic |

### Epigenomics - Spatial (`epigenomics`)

| Slug | Label | Resolution | Vendor |
|---|---|---|---|
| `spatial_atac` | Spatial ATAC-seq | single_cell | atlasxomics |
| `spatial_cut_tag` | Spatial CUT&Tag | single_cell | atlasxomics |
| `spatial_atac_rna` | Spatial ATAC-RNA-seq | single_cell | atlasxomics |
| `spatial_cut_tag_rna` | Spatial CUT&Tag-RNA-seq | single_cell | atlasxomics |
| `space_seq` | SPACE-seq | single_cell | atlasxomics |

### Genome Sequencing - Short Read (`genome_seq`)

| Slug | Label | Resolution | Vendor |
|---|---|---|---|
| `illumina_wgs` | Illumina WGS | bulk | illumina |
| `illumina_wes` | Illumina WES | bulk | illumina |
| `bgi_wgs` | BGI DNBSEQ WGS | bulk | bgi |
| `ultima_wgs` | Ultima Genomics UG 100 WGS | bulk | ultima |

### Genome Sequencing - Long Read (`genome_seq`)

| Slug | Label | Resolution | Vendor |
|---|---|---|---|
| `pacbio_hifi` | PacBio HiFi | bulk | pacbio |
| `ont_wgs` | Oxford Nanopore WGS | bulk | ont |
| `element_aviti` | Element Biosciences AVITI | bulk | element |

### Genome Sequencing - Targeted (`genome_seq`)

| Slug | Label | Resolution | Vendor |
|---|---|---|---|
| `targeted_panel` | Custom Gene Panels | bulk | generic |
| `amplicon_seq` | Amplicon Sequencing | bulk | generic |

---

## Usage Examples

### Look up a single platform

```python
from sc_tools.biodata import get_platform

spec = get_platform("xenium")
# PlatformSpec(name="xenium", biodata_type="spatial_seq",
#   modality="Spatial Transcriptomics - ISH", vendor="10x_genomics", ...)

spec.biodata_type   # "spatial_seq"
spec.modality       # "Spatial Transcriptomics - ISH"
spec.resolution     # "subcellular"
spec.defaults       # {"spatial_resolution": "single_cell", "coordinate_system": "micron"}
```

### List platforms with filters

```python
from sc_tools.biodata import list_platforms

# All spatial sequencing platforms
spatial_seq = list_platforms(category="spatial_seq")

# All 10x Genomics platforms
tenx = list_platforms(vendor="10x_genomics")

# All single-cell resolution platforms
sc = list_platforms(resolution="single_cell")

# Combine filters: spatial + protein measurement
spatial_prot = list_platforms(spatial=True, measurement="protein")
```

### List platforms by modality

```python
from sc_tools.biodata import list_platforms_by_modality

ish_platforms = list_platforms_by_modality("Spatial Transcriptomics - ISH")
# [PlatformSpec(name="cosmx_1k", ...), PlatformSpec(name="xenium", ...), ...]
```

### List available modalities

```python
from sc_tools.biodata import list_modalities

# All modalities
all_mods = list_modalities()
# ["Bulk RNA-seq", "Epigenomics - Bulk", "Epigenomics - Single Cell", ...]

# Modalities for a specific biodata_type
image_mods = list_modalities(biodata_type="image")
# ["Histology", "Spatial Proteomics - Cyclic IF",
#  "Spatial Proteomics - Mass Spec", "Spatial Proteomics - Standard IF/IHC"]
```

---

## Registering Custom Platforms

Use `register_platform()` to add project-specific or emerging platforms at runtime. No schema migration is needed -- the platform maps to an existing `BioDataType` ORM subclass.

```python
from sc_tools.biodata import PlatformSpec, register_platform

register_platform(
    "custom_panel",
    PlatformSpec(
        name="custom_panel",
        label="In-House Custom Panel v2",
        biodata_type="spatial_seq",
        category="spatial_seq",
        subcategory="imaging_based",
        measurement="rna",
        resolution="single_cell",
        spatial=True,
        vendor="in_house",
        modality="Spatial Transcriptomics - ISH",
        platform_version="v2",
        defaults={"panel_size": 500},
    ),
)
```

Key rules for custom platforms:

- `biodata_type` must be one of the five JTI types: `spatial_seq`, `image`, `rnaseq`, `epigenomics`, `genome_seq`.
- `modality` should match an existing modality string from the taxonomy table above, or introduce a new one if the platform represents a genuinely new experimental family.
- Calling `register_platform()` with an existing slug overwrites the previous entry.
- Custom registrations persist only for the current Python session. To make them permanent, add them to `sc_tools/biodata.py`.
