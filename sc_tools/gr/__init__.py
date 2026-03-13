"""
sc_tools.gr — multi-ROI wrappers around squidpy.gr functions.

All functions operate per-ROI and aggregate results across ROIs, ensuring
that spatial graphs are never computed across concatenated ROI boundaries.

Public API
----------
spatial_neighbors : Build block-diagonal spatial connectivity graph.
nhood_enrichment : Neighborhood enrichment per ROI with BH-corrected aggregation.
interaction_matrix : Cell-type interaction counts per ROI.
co_occurrence : Co-occurrence scores per ROI.
centrality_scores : Graph centrality scores per ROI.
ligrec : Ligand-receptor analysis per ROI.
ripley : Ripley statistics per ROI.
spatial_autocorr : Moran/Geary autocorrelation per ROI.
unify_matrices : Utility — align per-ROI matrices to global category space.
combine_pvalues : Utility — combine p-values across ROIs.
cross_roi_zscore : Utility — z-score values across ROIs.
"""

from ._aggregate import combine_pvalues, cross_roi_zscore, unify_matrices
from ._autocorr import spatial_autocorr
from ._centrality import centrality_scores
from ._graph import spatial_neighbors
from ._ligrec import ligrec
from ._nhood import interaction_matrix, nhood_enrichment
from ._occurrence import co_occurrence
from ._ripley import ripley

__all__ = [
    "spatial_neighbors",
    "nhood_enrichment",
    "interaction_matrix",
    "co_occurrence",
    "centrality_scores",
    "ligrec",
    "ripley",
    "spatial_autocorr",
    "unify_matrices",
    "combine_pvalues",
    "cross_roi_zscore",
]
