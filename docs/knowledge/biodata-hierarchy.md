# BioData Hierarchy

Full taxonomy of `BioDataType` → `BioDataModality` → `BioDataPlatform` used by `sc_tools.registry`.
For schema details, see [[registry]].

---

## Type Hierarchy (Joined Table Inheritance)

```
BioData (bio_data) ── base table
├── BioImage             ── imaging data without expression
├── RNASeqData           ── RNA expression (bulk or single-cell)
├── SpatialSeqData       ── spatially-resolved expression
├── EpigenomicsData      ── chromatin / epigenetic assays
└── GenomeSeqData        ── DNA sequencing
```

---

## Modalities and Platforms

### BioImage

| Modality | Platforms |
|----------|-----------|
| Imaging Mass Cytometry | `imc` (Fluidigm/Standard BioTools Hyperion) |
| Multiplexed Ion Beam Imaging | `mibi` |
| Mass Spectrometry Imaging | `maldi_ims` |
| Multiplexed IF | `codex`, `mibi_if`, `phenocycler` |
| Standard IF / IHC | `if_standard`, `ihc_standard` |
| Histology | `he_standard`, `he_wsi` |

### RNASeqData

| Modality | Platforms |
|----------|-----------|
| Single-cell RNA-seq | `10x_chromium_3p`, `10x_chromium_5p`, `dropseq`, `smartseq2`, `smartseq3` |
| Single-nucleus RNA-seq | `10x_chromium_snrna`, `parse_split_seq` |
| Bulk RNA-seq | `bulk_rnaseq_illumina`, `bulk_rnaseq_ont` |

### SpatialSeqData

| Modality | Platforms |
|----------|-----------|
| Visium (spot-level) | `visium` |
| Visium HD (bin-level) | `visium_hd` |
| Visium HD Cell (cell-segmented) | `visium_hd_cell` |
| Xenium | `xenium`, `xenium_5k` |
| CosMx | `cosmx_1k`, `cosmx_6k`, `cosmx_full_library` |
| MERFISH | `merfish` |
| seqFISH | `seqfish`, `seqfish_plus` |
| Slide-seq | `slideseq`, `slideseqv2` |
| Stereo-seq | `stereoseq` |
| HDST | `hdst` |

### EpigenomicsData

| Modality | Platforms |
|----------|-----------|
| ATAC-seq | `10x_atac`, `10x_multiome_atac`, `bulk_atac` |
| ChIP-seq | `chip_seq_illumina` |
| CUT&Tag / CUT&RUN | `cut_and_tag`, `cut_and_run` |
| Bisulfite / Methylation | `wgbs`, `rrbs` |

### GenomeSeqData

| Modality | Platforms |
|----------|-----------|
| Whole Genome Sequencing | `wgs_illumina`, `wgs_ont`, `wgs_pacbio` |
| Whole Exome Sequencing | `wes_illumina` |
| Targeted Panel | `targeted_panel_illumina` |

---

## Usage

```python
from sc_tools.registry import Registry

reg = Registry()

# Register a Visium HD dataset
reg.register_biodata(
    project_name="ggo_visium",
    sample_id="sample_001",
    uri="sftp://brb//athena/.../adata.filtered.h5ad",
    platform="visium_hd",       # auto-fills modality, subtype, bin_size
    phase="qc_filter",
)

# List all spatial data for a project
reg.list_biodata(project_name="ggo_visium", modality_family="Spatial Transcriptomics")
```

---

## Notes

- Modality family groups related platforms for cross-platform queries (e.g. all spatial transcriptomics).
- Platform slug is the lowest-level identifier; it uniquely determines the BioData subtype.
- `register_platform()` can extend the registry at runtime without schema migrations.
- Legacy `register_dataset()` dual-writes to `bio_data`; prefer `register_biodata()` for new data.
