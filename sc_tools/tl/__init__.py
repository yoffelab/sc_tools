"""
sc_tools.tl: Analysis tools for spatial omics data.

Following scanpy's API pattern, this module provides analysis functions:
- testing: Statistical testing (Mann-Whitney, FDR correction)
- colocalization: Spatial colocalization analysis (correlation, Moran's I, enrichment)
- deconvolution: Cell-type deconvolution (cell2location, tangram, destvi)
- gene_sets: Gene set loaders and curation utilities
- gsea: Group-level enrichment testing (ORA, pseudobulk GSEA)
"""

# Import tool functions
from . import colocalization, gene_sets, testing
from . import io as io_save
from .deconvolution import (
    deconvolution,
    extract_reference_profiles,
    select_signature_genes,
)
from .gene_sets import (
    list_gene_sets,
    load_gmt,
    load_hallmark,
    load_msigdb_json,
    merge_gene_signatures,
    save_gene_signatures,
    update_gene_symbols,
    validate_gene_signatures,
)
from .gsea import run_gsea_pseudobulk, run_ora
from .score_signature import score_signature

__all__ = [
    "testing",
    "colocalization",
    "io_save",
    "score_signature",
    "deconvolution",
    "extract_reference_profiles",
    "select_signature_genes",
    "gene_sets",
    "load_hallmark",
    "load_msigdb_json",
    "load_gmt",
    "list_gene_sets",
    "validate_gene_signatures",
    "merge_gene_signatures",
    "update_gene_symbols",
    "save_gene_signatures",
    "run_ora",
    "run_gsea_pseudobulk",
]
