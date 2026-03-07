"""
Shared figure configuration for DLBCL IMC manuscript.

Provides consistent colors, ordering, and matplotlib rcParams across all
figure scripts. Import this module in every fig*.py script:

    from figure_config import (
        LME_COLORS, LME_ORDER, CELLTYPE_COLORS, COO_COLORS,
        apply_figure_style,
    )
"""

import matplotlib as mpl

# ---------------------------------------------------------------------------
# LME classes — 5 classes, manuscript order (Cold first as largest group)
# ---------------------------------------------------------------------------

LME_ORDER = [
    "Cold",
    "Stromal",
    "Cytotoxic",
    "T cell Regulated",
    "CD206 Enriched",
]

LME_COLORS = {
    "Cold": "#4575b4",            # blue
    "Stromal": "#91bfdb",         # light blue
    "Cytotoxic": "#fc8d59",       # orange
    "T cell Regulated": "#fee090",  # yellow
    "CD206 Enriched": "#d73027",  # red
}

# ---------------------------------------------------------------------------
# Cell types — 12 major types from manuscript (immune + stromal panels)
# Okabe-Ito + Paul Tol extended for >8 categories
# ---------------------------------------------------------------------------

CELLTYPE_COLORS = {
    # Lymphoid
    "B cell": "#0072B2",
    "CD4 T cell": "#009E73",
    "CD8 T cell": "#56B4E9",
    "Treg": "#CC79A7",
    "NK": "#D55E00",
    # Myeloid
    "Macrophage": "#E69F00",
    "Monocyte": "#F0E442",
    "DC": "#999999",
    # Stromal
    "BEC": "#882255",
    "LEC": "#AA4499",
    "NESC": "#44AA99",
    # Other
    "Other": "#DDDDDD",
}

# Manuscript abbreviations for axis labels
CELLTYPE_DISPLAY = {
    "B cell": "B cell",
    "CD4 T cell": "CD4+ T",
    "CD8 T cell": "CD8+ T",
    "Treg": "Treg",
    "NK": "NK",
    "Macrophage": "Mac",
    "Monocyte": "Mono",
    "DC": "DC",
    "BEC": "BEC",
    "LEC": "LEC",
    "NESC": "NESC",
    "Other": "Other",
}

# ---------------------------------------------------------------------------
# Cell of Origin (COO)
# ---------------------------------------------------------------------------

COO_COLORS = {
    "GCB": "#4DAF4A",
    "ABC": "#E41A1C",
    "Unclassified": "#999999",
    "U": "#999999",  # alias in DLC380_clinical.tsv
}

# ---------------------------------------------------------------------------
# Matplotlib rcParams for publication figures
# ---------------------------------------------------------------------------

FIGURE_RC_PARAMS = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "legend.title_fontsize": 7,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "lines.linewidth": 0.8,
    "patch.linewidth": 0.5,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    "pdf.fonttype": 42,  # TrueType fonts in PDF (editable in Illustrator)
    "ps.fonttype": 42,
}


def apply_figure_style():
    """Apply publication-quality rcParams. Call once at script start."""
    mpl.rcParams.update(FIGURE_RC_PARAMS)


def get_lme_palette(classes=None):
    """Return list of colors matching LME_ORDER (or a subset)."""
    if classes is None:
        classes = LME_ORDER
    return [LME_COLORS.get(c, "#999999") for c in classes]


def significance_label(p: float) -> str:
    """Convert p-value to significance stars (BH-corrected input expected)."""
    if p < 0.0001:
        return "****"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


# ---------------------------------------------------------------------------
# Morphological / non-protein feature patterns (to filter from IMC var_names)
# ---------------------------------------------------------------------------

MORPHOLOGICAL_PATTERNS = [
    "solidity", "eccentricity", "filled-area", "perimeter", "euler",
    "extent", "area-", "major-axis", "minor-axis", "orientation",
    "centroid", "equivalent", "convex", "N-S", "full-index", "full_index",
    "nCount_Protein", "nFeature_Protein",
]


def is_protein_marker(name: str) -> bool:
    """Return True if marker name is a protein, not a morphological feature."""
    name_lower = name.lower()
    for pat in MORPHOLOGICAL_PATTERNS:
        if pat.lower() in name_lower:
            return False
    return True


def filter_protein_markers(var_names) -> list[str]:
    """Filter var_names to only protein markers (exclude morphological features)."""
    return [v for v in var_names if is_protein_marker(v)]


# ---------------------------------------------------------------------------
# Cell type label normalization (Seurat labels → display names)
# ---------------------------------------------------------------------------

# Map from common Seurat label patterns to CELLTYPE_COLORS keys
_CELLTYPE_ALIASES = {
    "bcell": "B cell",
    "b_cell": "B cell",
    "b cell": "B cell",
    "cd4_tcell": "CD4 T cell",
    "cd4 tcell": "CD4 T cell",
    "cd4 t cell": "CD4 T cell",
    "cd4_t_cell": "CD4 T cell",
    "cd4": "CD4 T cell",
    "cd8_tcell": "CD8 T cell",
    "cd8 tcell": "CD8 T cell",
    "cd8 t cell": "CD8 T cell",
    "cd8_t_cell": "CD8 T cell",
    "cd8": "CD8 T cell",
    "treg": "Treg",
    "t_reg": "Treg",
    "nk": "NK",
    "nk cell": "NK",
    "macrophage": "Macrophage",
    "mac": "Macrophage",
    "monocyte": "Monocyte",
    "mono": "Monocyte",
    "dc": "DC",
    "dendritic": "DC",
    "bec": "BEC",
    "endothelial": "BEC",
    "lec": "LEC",
    "nesc": "NESC",
    "stroma": "NESC",
    "stromal": "NESC",
    "fibroblast": "NESC",
    "tcell": "CD4 T cell",  # generic T cell → CD4 as default
    "other": "Other",
    "unknown": "Unknown",
}


def normalize_celltype(label: str) -> str:
    """Normalize a cell type label to manuscript display name."""
    label_str = str(label).strip()
    # Exact match in CELLTYPE_COLORS
    if label_str in CELLTYPE_COLORS:
        return label_str
    # Case-insensitive alias lookup
    label_lower = label_str.lower()
    if label_lower in _CELLTYPE_ALIASES:
        return _CELLTYPE_ALIASES[label_lower]
    # Partial match: check if any alias key is contained in label
    for alias, mapped in _CELLTYPE_ALIASES.items():
        if alias in label_lower:
            return mapped
    return label_str  # return as-is if no match


def build_celltype_palette(categories, fallback_cmap="tab20", normalize=True):
    """Build a color palette dict for cell type categories.

    Uses CELLTYPE_COLORS where possible, falls back to tab20 for unmapped.
    Set normalize=False to skip fuzzy matching (useful when you want distinct
    colors for fine-grained subtypes like T0_CD8+PD1+LAG3+).
    """
    import matplotlib.pyplot as plt

    palette = {}
    unmapped = []
    for cat in categories:
        cat_str = str(cat)
        # Try direct match
        if cat_str in CELLTYPE_COLORS:
            palette[cat_str] = CELLTYPE_COLORS[cat_str]
        elif normalize:
            # Try normalized match
            norm = normalize_celltype(cat_str)
            if norm in CELLTYPE_COLORS:
                palette[cat_str] = CELLTYPE_COLORS[norm]
            else:
                unmapped.append(cat_str)
        else:
            unmapped.append(cat_str)

    # Assign colors from combined tab20 + tab20b for unmapped categories
    if unmapped:
        cmap1 = plt.get_cmap("tab20")
        cmap2 = plt.get_cmap("tab20b")
        for i, cat in enumerate(unmapped):
            if i < 20:
                palette[cat] = cmap1(i / 20)
            else:
                palette[cat] = cmap2((i - 20) / 20)

    return palette


def is_bcell_label(label: str) -> bool:
    """Return True if label represents a B cell type."""
    label_lower = str(label).lower().strip()
    return any(pat in label_lower for pat in ["b cell", "b_cell", "bcell", "b-cell"])
