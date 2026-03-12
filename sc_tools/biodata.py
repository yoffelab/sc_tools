"""BioData platform registry for sc_tools.

Pure Python module (no DB dependency) providing a registry of known biological
data platforms. Each platform maps to a BioData ORM subtype with default field
values, enabling automatic classification when registering new data.

Usage::

    from sc_tools.biodata import get_platform, list_platforms, register_platform

    spec = get_platform("xenium")
    # PlatformSpec(name="xenium", biodata_type="spatial_seq", ...)

    spatial = list_platforms(category="spatial_seq")
    # [PlatformSpec(...), ...]
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass(frozen=True)
class PlatformSpec:
    """Specification for a known biological data platform."""

    name: str  # slug (e.g. "xenium")
    label: str  # human-readable (e.g. "10x Xenium")
    biodata_type: str  # "spatial_seq", "image", "rnaseq", "epigenomics", "genome_seq"
    category: str  # matches biodata_type
    subcategory: str  # e.g. "single_cell", "multiplexed", "spot"
    measurement: str  # "rna", "protein", "dna", "chromatin", "metabolite", "morphology"
    resolution: str  # "bulk", "spot", "single_cell", "subcellular"
    spatial: bool
    vendor: str  # "10x_genomics", "nanostring", etc.
    modality: str = ""  # e.g. "Spatial Proteomics - Mass Spec"
    platform_version: str = ""  # e.g. "v2", "v4"
    defaults: dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Platform registry (mutable at runtime via register_platform)
# ---------------------------------------------------------------------------

KNOWN_PLATFORMS: dict[str, PlatformSpec] = {}


def register_platform(name: str, spec: PlatformSpec) -> None:
    """Register or overwrite a platform specification."""
    KNOWN_PLATFORMS[name] = spec


def get_platform(name: str) -> PlatformSpec:
    """Return the PlatformSpec for a known platform slug.

    Raises KeyError if the platform is not registered.
    """
    if name not in KNOWN_PLATFORMS:
        raise KeyError(
            f"Unknown platform '{name}'. "
            f"Use register_platform() to add it or check list_platforms()."
        )
    return KNOWN_PLATFORMS[name]


def list_platforms(
    category: str | None = None,
    vendor: str | None = None,
    spatial: bool | None = None,
    measurement: str | None = None,
    resolution: str | None = None,
) -> list[PlatformSpec]:
    """Return platform specs matching the given filters.

    All filters are optional; when omitted that dimension is not filtered.
    """
    results = []
    for spec in KNOWN_PLATFORMS.values():
        if category is not None and spec.category != category:
            continue
        if vendor is not None and spec.vendor != vendor:
            continue
        if spatial is not None and spec.spatial != spatial:
            continue
        if measurement is not None and spec.measurement != measurement:
            continue
        if resolution is not None and spec.resolution != resolution:
            continue
        results.append(spec)
    return sorted(results, key=lambda s: s.name)


def platform_for_project(project_platform: str) -> PlatformSpec | None:
    """Map an existing project platform string to a PlatformSpec.

    Returns None if no match is found (caller should handle gracefully).
    """
    return KNOWN_PLATFORMS.get(project_platform)


def list_modalities(biodata_type: str | None = None) -> list[str]:
    """Return sorted unique modality strings, optionally filtered by biodata_type."""
    modalities: set[str] = set()
    for spec in KNOWN_PLATFORMS.values():
        if not spec.modality:
            continue
        if biodata_type is not None and spec.biodata_type != biodata_type:
            continue
        modalities.add(spec.modality)
    return sorted(modalities)


def list_platforms_by_modality(modality: str) -> list[PlatformSpec]:
    """Return all platforms belonging to the given modality string."""
    return sorted(
        (s for s in KNOWN_PLATFORMS.values() if s.modality == modality),
        key=lambda s: s.name,
    )


def get_modality_for_platform(name: str) -> str:
    """Return the modality string for a platform slug. Raises KeyError if unknown."""
    return get_platform(name).modality


# ---------------------------------------------------------------------------
# Pre-seed: SPATIAL TRANSCRIPTOMICS — Sequencing-based (spot resolution)
# ---------------------------------------------------------------------------

_SPATIAL_SEQ_SPOT = [
    ("visium", "10x Visium", "spot", "10x_genomics", {"spatial_resolution": "spot"}),
    (
        "visium_hd",
        "10x Visium HD",
        "spot",
        "10x_genomics",
        {"spatial_resolution": "spot", "bin_size_um": 8.0},
    ),
    ("visium_hd_cell", "10x Visium HD Cell", "single_cell", "10x_genomics", {}),
    ("stereo_seq", "BGI Stereo-seq", "subcellular", "bgi", {}),
    ("slide_seq", "Slide-seq / Slide-seqV2", "spot", "broad_institute", {}),
    ("hdst", "High-Definition Spatial Transcriptomics", "spot", "broad_institute", {}),
    ("pixel_seq", "PIXEL-seq", "subcellular", "academic", {}),
    ("seq_scope", "Seq-Scope", "subcellular", "academic", {}),
    ("open_st", "Open-ST", "subcellular", "academic", {}),
    ("dbit_seq", "DBiT-seq", "spot", "academic", {}),
    ("sci_space", "sci-Space", "spot", "academic", {}),
    ("magic_seq", "MAGIC-seq", "spot", "academic", {}),
]

for _name, _label, _res, _vendor, _defaults in _SPATIAL_SEQ_SPOT:
    register_platform(
        _name,
        PlatformSpec(
            name=_name,
            label=_label,
            biodata_type="spatial_seq",
            category="spatial_seq",
            subcategory="sequencing_based",
            measurement="rna",
            resolution=_res,
            spatial=True,
            vendor=_vendor,
            modality="Spatial Transcriptomics - Sequencing",
            defaults=_defaults,
        ),
    )

# ---------------------------------------------------------------------------
# Pre-seed: SPATIAL TRANSCRIPTOMICS — Imaging-based / ISH
# ---------------------------------------------------------------------------

_SPATIAL_SEQ_ISH = [
    (
        "xenium",
        "10x Xenium",
        "subcellular",
        "10x_genomics",
        {"spatial_resolution": "single_cell", "coordinate_system": "micron"},
    ),
    ("cosmx_1k", "NanoString CosMx SMI 1k", "single_cell", "nanostring", {"panel_size": 1000}),
    ("cosmx_6k", "NanoString CosMx SMI 6k", "single_cell", "nanostring", {"panel_size": 6000}),
    ("cosmx_wt", "NanoString CosMx Whole Transcriptome", "single_cell", "nanostring", {}),
    ("merscope", "Vizgen MERSCOPE / MERFISH", "subcellular", "vizgen", {}),
    ("seqfish", "seqFISH / seqFISH+", "subcellular", "academic", {}),
    ("starmap", "STARmap / STARmap PLUS", "subcellular", "academic", {}),
    ("hybiss", "HybISS", "subcellular", "academic", {}),
    ("eelfish", "EEL FISH", "subcellular", "academic", {}),
    ("resolve", "Resolve Biosciences Molecular Cartography", "subcellular", "resolve", {}),
]

for _name, _label, _res, _vendor, _defaults in _SPATIAL_SEQ_ISH:
    register_platform(
        _name,
        PlatformSpec(
            name=_name,
            label=_label,
            biodata_type="spatial_seq",
            category="spatial_seq",
            subcategory="imaging_based",
            measurement="rna",
            resolution=_res,
            spatial=True,
            vendor=_vendor,
            modality="Spatial Transcriptomics - ISH",
            defaults=_defaults,
        ),
    )

# ---------------------------------------------------------------------------
# Pre-seed: SPATIAL TRANSCRIPTOMICS — Region-capture
# ---------------------------------------------------------------------------

for _name, _label, _vendor in [
    ("geomx", "NanoString GeoMx DSP", "nanostring"),
    ("visium_cytassist", "10x CytAssist", "10x_genomics"),
]:
    register_platform(
        _name,
        PlatformSpec(
            name=_name,
            label=_label,
            biodata_type="spatial_seq",
            category="spatial_seq",
            subcategory="region_capture",
            measurement="rna",
            resolution="spot",
            spatial=True,
            vendor=_vendor,
            modality="Spatial Transcriptomics - Region Capture",
        ),
    )

# ---------------------------------------------------------------------------
# Pre-seed: SPATIAL PROTEOMICS — Mass spectrometry imaging
# ---------------------------------------------------------------------------

_SPATIAL_PROT_MS = [
    (
        "imc",
        "Imaging Mass Cytometry",
        "multiplexed",
        "standard_biotools",
        {"staining_protocol": "IMC", "image_type": "multiplexed"},
    ),
    (
        "mibi",
        "Multiplexed Ion Beam Imaging",
        "multiplexed",
        "ionpath",
        {"staining_protocol": "MIBI", "image_type": "multiplexed"},
    ),
    ("maldi_ims", "MALDI Imaging Mass Spec", "bulk", "bruker", {}),
]

for _name, _label, _subcat, _vendor, _defaults in _SPATIAL_PROT_MS:
    register_platform(
        _name,
        PlatformSpec(
            name=_name,
            label=_label,
            biodata_type="image",
            category="image",
            subcategory=_subcat,
            measurement="protein",
            resolution="single_cell",
            spatial=True,
            vendor=_vendor,
            modality="Spatial Proteomics - Mass Spec",
            defaults=_defaults,
        ),
    )

# ---------------------------------------------------------------------------
# Pre-seed: SPATIAL PROTEOMICS — Cyclic immunofluorescence
# ---------------------------------------------------------------------------

_SPATIAL_PROT_CIF = [
    (
        "phenocycler",
        "Akoya PhenoCycler (CODEX)",
        "akoya",
        {"staining_protocol": "CODEX", "image_type": "multiplexed"},
    ),
    (
        "macsima",
        "Miltenyi MACSima",
        "miltenyi",
        {"staining_protocol": "MACSima", "image_type": "multiplexed"},
    ),
    ("cycif", "CyCIF", "academic", {"staining_protocol": "CyCIF", "image_type": "multiplexed"}),
    ("insituplex", "Ultivue InSituPlex", "ultivue", {"image_type": "multiplexed"}),
    ("comet", "Lunaphore COMET", "lunaphore", {"image_type": "multiplexed"}),
    ("orion", "RareCyte ORION", "rarecyte", {"image_type": "multiplexed"}),
]

for _name, _label, _vendor, _defaults in _SPATIAL_PROT_CIF:
    register_platform(
        _name,
        PlatformSpec(
            name=_name,
            label=_label,
            biodata_type="image",
            category="image",
            subcategory="cyclic_if",
            measurement="protein",
            resolution="single_cell",
            spatial=True,
            vendor=_vendor,
            modality="Spatial Proteomics - Cyclic IF",
            defaults=_defaults,
        ),
    )

# ---------------------------------------------------------------------------
# Pre-seed: SPATIAL PROTEOMICS — Standard IF/IHC
# ---------------------------------------------------------------------------

for _name, _label, _vendor in [
    ("mihc", "Multiplex IHC (Vectra/Polaris)", "akoya"),
    ("if_standard", "Standard Immunofluorescence", "generic"),
    ("vectra_polaris", "Akoya Vectra Polaris", "akoya"),
]:
    register_platform(
        _name,
        PlatformSpec(
            name=_name,
            label=_label,
            biodata_type="image",
            category="image",
            subcategory="standard_if_ihc",
            measurement="protein",
            resolution="single_cell",
            spatial=True,
            vendor=_vendor,
            modality="Spatial Proteomics - Standard IF/IHC",
        ),
    )

# ---------------------------------------------------------------------------
# Pre-seed: HISTOLOGY / PATHOLOGY IMAGES
# ---------------------------------------------------------------------------

_HISTOLOGY = [
    ("he", "H&E Staining", "he"),
    ("brightfield", "Unstained Brightfield", "brightfield"),
    ("phase_contrast", "Phase Contrast Microscopy", "phase_contrast"),
    ("pas", "PAS Staining", "pas"),
    ("masson_trichrome", "Masson Trichrome", "masson_trichrome"),
    ("electron_microscopy", "Electron Microscopy (TEM/SEM)", "electron_microscopy"),
]

for _name, _label, _img_type in _HISTOLOGY:
    register_platform(
        _name,
        PlatformSpec(
            name=_name,
            label=_label,
            biodata_type="image",
            category="image",
            subcategory="histology",
            measurement="morphology",
            resolution="bulk",
            spatial=True,
            vendor="generic",
            modality="Histology",
            defaults={
                "image_type": _img_type,
                "staining_protocol": "H&E" if _name == "he" else _img_type,
            },
        ),
    )

# ---------------------------------------------------------------------------
# Pre-seed: SINGLE-CELL RNA-seq — Droplet-based
# ---------------------------------------------------------------------------

_SCRNA_DROPLET = [
    ("chromium_3p", "10x Chromium 3'", "chromium_v3", "10x_genomics"),
    ("chromium_5p", "10x Chromium 5'", "chromium_v3", "10x_genomics"),
    ("chromium_flex", "10x Chromium Flex", "chromium_flex", "10x_genomics"),
    ("chromium_multiome", "10x Multiome (RNA + ATAC)", "chromium_multiome", "10x_genomics"),
    ("indrops", "inDrops", "indrops", "academic"),
    ("dropseq", "Drop-seq", "dropseq", "academic"),
]

for _name, _label, _chem, _vendor in _SCRNA_DROPLET:
    register_platform(
        _name,
        PlatformSpec(
            name=_name,
            label=_label,
            biodata_type="rnaseq",
            category="rnaseq",
            subcategory="droplet",
            measurement="rna",
            resolution="single_cell",
            spatial=False,
            vendor=_vendor,
            modality="Single-Cell RNA-seq - Droplet",
            defaults={"chemistry": _chem, "library_type": "single_cell"},
        ),
    )

# ---------------------------------------------------------------------------
# Pre-seed: SINGLE-CELL RNA-seq — Combinatorial barcoding
# ---------------------------------------------------------------------------

for _name, _label, _chem, _vendor in [
    ("parse_evercode", "Parse Biosciences Evercode", "evercode_v3", "parse"),
    ("scale_rna", "Scale Biosciences", "scale", "scale"),
    ("split_seq", "SPLiT-seq", "split_seq", "academic"),
]:
    register_platform(
        _name,
        PlatformSpec(
            name=_name,
            label=_label,
            biodata_type="rnaseq",
            category="rnaseq",
            subcategory="combinatorial",
            measurement="rna",
            resolution="single_cell",
            spatial=False,
            vendor=_vendor,
            modality="Single-Cell RNA-seq - Combinatorial",
            defaults={"chemistry": _chem, "library_type": "single_cell"},
        ),
    )

# ---------------------------------------------------------------------------
# Pre-seed: SINGLE-CELL RNA-seq — Microwell
# ---------------------------------------------------------------------------

for _name, _label, _chem, _vendor in [
    ("bd_rhapsody", "BD Rhapsody", "rhapsody", "bd"),
    ("seq_well", "Seq-Well", "seq_well", "academic"),
    ("microwell_seq", "Microwell-seq", "microwell", "academic"),
]:
    register_platform(
        _name,
        PlatformSpec(
            name=_name,
            label=_label,
            biodata_type="rnaseq",
            category="rnaseq",
            subcategory="microwell",
            measurement="rna",
            resolution="single_cell",
            spatial=False,
            vendor=_vendor,
            modality="Single-Cell RNA-seq - Microwell",
            defaults={"chemistry": _chem, "library_type": "single_cell"},
        ),
    )

# ---------------------------------------------------------------------------
# Pre-seed: SINGLE-CELL RNA-seq — Plate-based
# ---------------------------------------------------------------------------

for _name, _label, _chem in [
    ("smart_seq2", "Smart-seq2", "smart_seq2"),
    ("smart_seq3", "Smart-seq3", "smart_seq3"),
    ("cel_seq2", "CEL-Seq2", "cel_seq2"),
]:
    register_platform(
        _name,
        PlatformSpec(
            name=_name,
            label=_label,
            biodata_type="rnaseq",
            category="rnaseq",
            subcategory="plate_based",
            measurement="rna",
            resolution="single_cell",
            spatial=False,
            vendor="academic",
            modality="Single-Cell RNA-seq - Plate",
            defaults={"chemistry": _chem, "library_type": "single_cell"},
        ),
    )

# ---------------------------------------------------------------------------
# Pre-seed: BULK RNA-seq
# ---------------------------------------------------------------------------

for _name, _label, _vendor, _defaults in [
    (
        "illumina_rnaseq",
        "Illumina Short-Read RNA-seq",
        "illumina",
        {"sequencing_platform": "illumina_novaseq", "library_type": "bulk"},
    ),
    (
        "bgi_rnaseq",
        "BGI DNBSEQ RNA-seq",
        "bgi",
        {"sequencing_platform": "bgi_dnbseq", "library_type": "bulk"},
    ),
    (
        "ont_direct_rna",
        "Oxford Nanopore Direct RNA-seq",
        "ont",
        {"sequencing_platform": "ont_promethion", "library_type": "bulk"},
    ),
    (
        "pacbio_isoseq",
        "PacBio Iso-Seq",
        "pacbio",
        {"sequencing_platform": "pacbio_hifi", "library_type": "bulk"},
    ),
    (
        "ultima_rnaseq",
        "Ultima Genomics UG 100",
        "ultima",
        {"sequencing_platform": "ultima_ug100", "library_type": "bulk"},
    ),
]:
    register_platform(
        _name,
        PlatformSpec(
            name=_name,
            label=_label,
            biodata_type="rnaseq",
            category="rnaseq",
            subcategory="bulk",
            measurement="rna",
            resolution="bulk",
            spatial=False,
            vendor=_vendor,
            modality="Bulk RNA-seq",
            defaults=_defaults,
        ),
    )

# ---------------------------------------------------------------------------
# Pre-seed: SINGLE-CELL MULTIOMICS
# ---------------------------------------------------------------------------

_MULTIOMICS = [
    ("cite_seq", "CITE-seq", "rna", "biolegend"),
    ("reap_seq", "REAP-seq", "rna", "academic"),
    ("asap_seq", "ASAP-seq", "chromatin", "academic"),
    ("share_seq", "SHARE-seq", "rna", "academic"),
    ("dogma_seq", "DOGMA-seq", "rna", "academic"),
    ("tea_seq", "TEA-seq", "rna", "academic"),
    ("snare_seq", "SNARE-seq", "rna", "academic"),
    ("paired_seq", "Paired-seq", "rna", "academic"),
]

for _name, _label, _meas, _vendor in _MULTIOMICS:
    register_platform(
        _name,
        PlatformSpec(
            name=_name,
            label=_label,
            biodata_type="rnaseq",
            category="rnaseq",
            subcategory="multiomics",
            measurement=_meas,
            resolution="single_cell",
            spatial=False,
            vendor=_vendor,
            modality="Multiomics",
        ),
    )

# ---------------------------------------------------------------------------
# Pre-seed: EPIGENOMICS — Bulk
# ---------------------------------------------------------------------------

_EPIGEN_BULK = [
    ("atac_seq", "ATAC-seq", "atac_seq", None),
    ("chip_seq", "ChIP-seq", "chip_seq", None),
    ("cut_and_run", "CUT&RUN", "cut_and_run", None),
    ("cut_and_tag", "CUT&Tag", "cut_and_tag", None),
    ("bisulfite_seq", "Bisulfite Sequencing (WGBS/RRBS)", "methylation", None),
    ("em_seq", "EM-seq", "methylation", None),
    ("hi_c", "Hi-C", "hi_c", None),
]

for _name, _label, _assay, _target in _EPIGEN_BULK:
    register_platform(
        _name,
        PlatformSpec(
            name=_name,
            label=_label,
            biodata_type="epigenomics",
            category="epigenomics",
            subcategory="bulk",
            measurement="chromatin",
            resolution="bulk",
            spatial=False,
            vendor="generic",
            modality="Epigenomics - Bulk",
            defaults={"assay_type": _assay},
        ),
    )

# ---------------------------------------------------------------------------
# Pre-seed: EPIGENOMICS — Single-cell
# ---------------------------------------------------------------------------

for _name, _label, _assay in [
    ("sc_atac_seq", "scATAC-seq", "atac_seq"),
    ("sc_cut_and_tag", "Single-cell CUT&Tag", "cut_and_tag"),
    ("sc_methylation", "Single-cell Bisulfite (scBS-seq)", "methylation"),
    ("sc_hi_c", "Single-cell Hi-C", "hi_c"),
]:
    register_platform(
        _name,
        PlatformSpec(
            name=_name,
            label=_label,
            biodata_type="epigenomics",
            category="epigenomics",
            subcategory="single_cell",
            measurement="chromatin",
            resolution="single_cell",
            spatial=False,
            vendor="generic",
            modality="Epigenomics - Single Cell",
            defaults={"assay_type": _assay},
        ),
    )

# ---------------------------------------------------------------------------
# Pre-seed: EPIGENOMICS — Spatial
# ---------------------------------------------------------------------------

for _name, _label, _assay in [
    ("spatial_atac", "Spatial ATAC-seq", "spatial_atac"),
    ("spatial_cut_tag", "Spatial CUT&Tag", "spatial_cut_tag"),
    ("spatial_atac_rna", "Spatial ATAC-RNA-seq", "spatial_atac"),
    ("spatial_cut_tag_rna", "Spatial CUT&Tag-RNA-seq", "spatial_cut_tag"),
    ("space_seq", "SPACE-seq", "space_seq"),
]:
    register_platform(
        _name,
        PlatformSpec(
            name=_name,
            label=_label,
            biodata_type="epigenomics",
            category="epigenomics",
            subcategory="spatial",
            measurement="chromatin",
            resolution="single_cell",
            spatial=True,
            vendor="atlasxomics",
            modality="Epigenomics - Spatial",
            defaults={"assay_type": _assay},
        ),
    )

# ---------------------------------------------------------------------------
# Pre-seed: GENOME SEQUENCING — Short-read
# ---------------------------------------------------------------------------

for _name, _label, _seq_type, _vendor in [
    ("illumina_wgs", "Illumina WGS", "wgs", "illumina"),
    ("illumina_wes", "Illumina WES", "wes", "illumina"),
    ("bgi_wgs", "BGI DNBSEQ WGS", "wgs", "bgi"),
    ("ultima_wgs", "Ultima Genomics UG 100 WGS", "wgs", "ultima"),
]:
    register_platform(
        _name,
        PlatformSpec(
            name=_name,
            label=_label,
            biodata_type="genome_seq",
            category="genome_seq",
            subcategory="short_read",
            measurement="dna",
            resolution="bulk",
            spatial=False,
            vendor=_vendor,
            modality="Genome Sequencing - Short Read",
            defaults={
                "sequencing_type": _seq_type,
                "sequencing_platform": f"{_vendor}_novaseq"
                if _vendor == "illumina"
                else f"{_vendor}_dnbseq"
                if _vendor == "bgi"
                else f"{_vendor}_ug100",
            },
        ),
    )

# ---------------------------------------------------------------------------
# Pre-seed: GENOME SEQUENCING — Long-read
# ---------------------------------------------------------------------------

for _name, _label, _vendor in [
    ("pacbio_hifi", "PacBio HiFi", "pacbio"),
    ("ont_wgs", "Oxford Nanopore WGS", "ont"),
    ("element_aviti", "Element Biosciences AVITI", "element"),
]:
    register_platform(
        _name,
        PlatformSpec(
            name=_name,
            label=_label,
            biodata_type="genome_seq",
            category="genome_seq",
            subcategory="long_read",
            measurement="dna",
            resolution="bulk",
            spatial=False,
            vendor=_vendor,
            modality="Genome Sequencing - Long Read",
            defaults={
                "sequencing_type": "wgs",
                "sequencing_platform": f"{_vendor}_hifi"
                if _vendor == "pacbio"
                else f"{_vendor}_promethion"
                if _vendor == "ont"
                else f"{_vendor}_aviti",
            },
        ),
    )

# ---------------------------------------------------------------------------
# Pre-seed: GENOME SEQUENCING — Targeted
# ---------------------------------------------------------------------------

for _name, _label in [
    ("targeted_panel", "Custom Gene Panels"),
    ("amplicon_seq", "Amplicon Sequencing"),
]:
    register_platform(
        _name,
        PlatformSpec(
            name=_name,
            label=_label,
            biodata_type="genome_seq",
            category="genome_seq",
            subcategory="targeted",
            measurement="dna",
            resolution="bulk",
            spatial=False,
            vendor="generic",
            modality="Genome Sequencing - Targeted",
            defaults={"sequencing_type": "targeted_panel"},
        ),
    )


# ---------------------------------------------------------------------------
# Clean up module-level loop variables
# ---------------------------------------------------------------------------
del _name, _label
# Other loop vars already cleaned by scope

__all__ = [
    "PlatformSpec",
    "KNOWN_PLATFORMS",
    "register_platform",
    "get_platform",
    "list_platforms",
    "platform_for_project",
    "list_modalities",
    "list_platforms_by_modality",
    "get_modality_for_platform",
]
