"""
Function to create a TLS-clustered anndata object from spatial Visium data.

Each observation in the output represents a spatially connected TLS structure,
with aggregated expression values.
"""

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as sp
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import connected_components
from typing import Optional, Union
import warnings

# Visium spot size: 55 μm diameter = ~0.002375 mm² per spot
VISIUM_SPOT_AREA_MM2 = 0.002375


def create_tls_anndata(
    adata: ad.AnnData,
    tls_col: str = "architecture_type",
    tls_pattern: str = "TLS",
    library_id_key: str = "library_id",
    spatial_key: str = "spatial",
    solidity_type_col: str = "solidity_type",
    connectivity_threshold: Optional[float] = None,
    use_raw: bool = True,
    n_neighbors: int = 6,
) -> ad.AnnData:
    """
    Create an anndata object where each observation is a spatially connected TLS cluster.
    
    Parameters
    ----------
    adata : AnnData
        Input anndata object with spatial coordinates and TLS annotations
    tls_col : str
        Column name in adata.obs that contains TLS annotations (default: "architecture_type")
    tls_pattern : str
        Pattern to match TLS annotations. Spots where tls_col contains this pattern
        will be considered TLS (default: "TLS")
    library_id_key : str
        Column name for library/sample ID (default: "library_id")
    spatial_key : str
        Key in adata.obsm containing spatial coordinates (default: "spatial")
    solidity_type_col : str
        Column name for solidity_type information (default: "solidity_type")
    connectivity_threshold : float, optional
        Maximum distance for two spots to be considered connected (in coordinate units).
        If None, uses n_neighbors to determine connectivity (default: None)
    use_raw : bool
        If True, uses adata.raw for expression data. Otherwise uses adata.X (default: True)
    n_neighbors : int
        Number of nearest neighbors to consider for connectivity if connectivity_threshold
        is None (default: 6)
    
    Returns
    -------
    AnnData
        New anndata object where:
        - Each observation is a unique TLS cluster
        - obs contains: library_id, tls_id (unique identifier), tls_size (number of spots),
          tls_size_mm2 (approximate size in mm²), solidity_type_majority, solidity_type_proportion
        - X contains mean gene expression across spots in each TLS
        - obsm['summed_expression'] contains summed expression (raw counts) per TLS
    """
    
    # Validate inputs
    assert tls_col in adata.obs.columns, f"Column '{tls_col}' not found in adata.obs"
    assert library_id_key in adata.obs.columns, f"Column '{library_id_key}' not found in adata.obs"
    assert spatial_key in adata.obsm, f"'{spatial_key}' not found in adata.obsm"
    
    # Check if solidity_type column exists (optional, but warn if missing)
    has_solidity = solidity_type_col in adata.obs.columns

    if not has_solidity:
        warnings.warn(f"Column '{solidity_type_col}' not found in adata.obs. Solidity information will not be included.")
    
    if use_raw:
        assert adata.raw is not None, "adata.raw is None but use_raw=True"
        X_expr = adata.raw.X
        var_names = adata.raw.var_names
    else:
        X_expr = adata.X
        var_names = adata.var_names
    
    # Filter TLS spots
    is_tls = adata.obs[tls_col].astype(str).str.contains(tls_pattern, case=False, na=False).values
    if is_tls.sum() == 0:
        warnings.warn(f"No TLS spots found with pattern '{tls_pattern}' in column '{tls_col}'")
        # Return empty anndata with same structure
        empty_obs = pd.DataFrame(columns=[library_id_key, 'tls_id', 'tls_size'])
        return ad.AnnData(obs=empty_obs, var=pd.DataFrame(index=var_names))
    
    adata_tls = adata[is_tls].copy()
    
    # Get spatial coordinates
    spatial_coords = adata_tls.obsm[spatial_key]
    library_ids = adata_tls.obs[library_id_key].values
    
    # Get solidity_type information if available
    if has_solidity:
        solidity_types = adata_tls.obs[solidity_type_col].values
        # Define tie-breaking order: Solid > Non-Solid > Normal (and any others)
        # This ensures consistent ordering when there are ties
        solidity_tie_order = ["Solid", "Non-Solid", "Normal"]
    
    # Initialize TLS cluster labels
    n_spots = adata_tls.n_obs
    tls_labels = np.full(n_spots, -1, dtype=int)
    
    # Process each library separately to avoid cross-sample connections
    tls_counter = 0
    
    for lib_id in np.unique(library_ids):
        lib_mask = library_ids == lib_id
        lib_indices = np.where(lib_mask)[0]
        
        if len(lib_indices) == 0:
            continue
        
        lib_coords = spatial_coords[lib_indices]
        
        # Handle single spot case
        if len(lib_indices) == 1:
            tls_labels[lib_indices[0]] = tls_counter
            tls_counter += 1
            continue
        
        # Build connectivity graph
        if connectivity_threshold is not None:
            # Use distance threshold
            tree = cKDTree(lib_coords)
            pairs = tree.query_pairs(connectivity_threshold, output_type='ndarray')
            
            # Build sparse adjacency matrix
            n_lib = len(lib_indices)
            adj_matrix = sp.lil_matrix((n_lib, n_lib), dtype=bool)
            if len(pairs) > 0:
                adj_matrix[pairs[:, 0], pairs[:, 1]] = True
                adj_matrix[pairs[:, 1], pairs[:, 0]] = True  # Make symmetric
            
        else:
            # Use k-nearest neighbors
            tree = cKDTree(lib_coords)
            distances, neighbors = tree.query(lib_coords, k=min(n_neighbors + 1, len(lib_indices)))
            
            # Build sparse adjacency matrix
            n_lib = len(lib_indices)
            adj_matrix = sp.lil_matrix((n_lib, n_lib), dtype=bool)
            
            for i in range(n_lib):
                # Skip self (first neighbor is always self with distance 0)
                for j in neighbors[i, 1:]:
                    if j < n_lib:  # Valid index
                        adj_matrix[i, j] = True
                        adj_matrix[j, i] = True  # Make symmetric
        
        # Find connected components
        adj_matrix = adj_matrix.tocsr()
        n_components, component_labels = connected_components(
            csgraph=adj_matrix, directed=False, return_labels=True
        )
        
        # Assign unique TLS IDs across all libraries
        for comp_id in range(n_components):
            comp_mask = component_labels == comp_id
            global_indices = lib_indices[comp_mask]
            tls_labels[global_indices] = tls_counter
            tls_counter += 1
    
    # Create output data structures
    unique_tls_ids = np.unique(tls_labels)
    n_tls = len(unique_tls_ids)
    
    # Initialize output arrays
    obs_list = []
    X_mean_list = []
    X_sum_list = []
    spatial_centroids = []
    
    # Get expression data for TLS spots (already filtered, so indices align)
    X_tls = X_expr[is_tls]
    
    for tls_id in unique_tls_ids:
        tls_mask = tls_labels == tls_id
        tls_spots = np.where(tls_mask)[0]
        
        # Get library_id (should be same for all spots in cluster, but take first)
        lib_id = library_ids[tls_spots[0]]
        
        # Compute spatial centroid (mean coordinates of all spots in this TLS)
        tls_coords = spatial_coords[tls_mask]
        centroid = np.mean(tls_coords, axis=0)
        spatial_centroids.append(centroid)
        
        # Compute solidity_type majority and proportion
        if has_solidity:
            tls_solidity = solidity_types[tls_mask]
            # Count occurrences of each solidity type
            unique, counts = np.unique(tls_solidity, return_counts=True)
            solidity_counts = dict(zip(unique, counts))
            
            # Find majority type with tie-breaking
            # Sort by count (descending), then by tie-breaking order
            max_count = max(counts)
            candidates = [sol for sol, count in solidity_counts.items() if count == max_count]
            
            # If tie, use tie-breaking order
            if len(candidates) > 1:
                # Sort candidates by tie-breaking order (earlier in list = higher priority)
                candidates_ordered = sorted(
                    candidates,
                    key=lambda x: (
                        solidity_tie_order.index(x) if x in solidity_tie_order else len(solidity_tie_order)
                    )
                )
                majority_solidity = candidates_ordered[0]
            else:
                majority_solidity = candidates[0]
            
            # Calculate proportion
            total_spots = len(tls_spots)
            majority_count = solidity_counts[majority_solidity]
            solidity_proportion = majority_count / total_spots
        else:
            majority_solidity = None
            solidity_proportion = None
        
        # Aggregate expression (tls_mask already aligns with X_tls since both are filtered)
        X_tls_subset = X_tls[tls_mask]
        
        # Mean expression (for X)
        if sp.issparse(X_tls_subset):
            X_mean = np.asarray(X_tls_subset.mean(axis=0)).ravel()
            X_sum = np.asarray(X_tls_subset.sum(axis=0)).ravel()
        else:
            X_mean = X_tls_subset.mean(axis=0)
            X_sum = X_tls_subset.sum(axis=0)
        
        # Calculate TLS size in mm²
        n_spots = len(tls_spots)
        tls_size_mm2 = n_spots * VISIUM_SPOT_AREA_MM2
        
        # Store results
        obs_dict = {
            library_id_key: lib_id,
            'tls_id': f"{lib_id}_TLS_{tls_id:04d}",
            'tls_size': n_spots,  # number of spots
            'tls_size_mm2': tls_size_mm2,  # approximate size in mm²
        }
        
        if has_solidity:
            obs_dict['solidity_type_majority'] = majority_solidity
            obs_dict['solidity_type_proportion'] = solidity_proportion
        
        obs_list.append(obs_dict)
        X_mean_list.append(X_mean)
        X_sum_list.append(X_sum)
    
    # Create output anndata
    obs_df = pd.DataFrame(obs_list)
    obs_df.index = obs_df['tls_id'].values
    
    X_mean_array = np.vstack(X_mean_list)
    X_sum_array = np.vstack(X_sum_list)
    
    # Create var dataframe (copy metadata if available)
    if use_raw and adata.raw is not None:
        var_df = adata.raw.var.copy()
    else:
        var_df = adata.var.copy()
    
    # Ensure index matches var_names
    var_df = var_df.reindex(var_names)
    
    # Create new anndata
    adata_tls_out = ad.AnnData(
        X=X_mean_array,
        obs=obs_df,
        var=var_df,
        dtype=X_mean_array.dtype
    )
    
    # Add summed expression to obsm
    adata_tls_out.obsm['summed_expression'] = X_sum_array
    
    # Add spatial centroids to obsm
    adata_tls_out.obsm[spatial_key] = np.vstack(spatial_centroids)
    
    # Copy over uns metadata if needed
    if 'spatial' in adata.uns:
        adata_tls_out.uns['spatial'] = adata.uns['spatial']
    
    return adata_tls_out


# Example usage:
if __name__ == "__main__":
    # Load data
    adata = sc.read('results/adata.normalized.scored.p35.h5ad')
    
    # Create TLS anndata
    adata_tls = create_tls_anndata(
        adata,
        tls_col="architecture_type",
        tls_pattern="TLS",
        connectivity_threshold=None,  # Use k-nearest neighbors (default: 6)
        use_raw=True,
    )
    
    print(f"Input: {adata.n_obs} spots")
    print(f"Output: {adata_tls.n_obs} TLS clusters")
    print(f"\nTLS clusters per library:")
    print(adata_tls.obs.groupby('library_id').size())
    print(f"\nTLS size distribution:")
    print(adata_tls.obs['tls_size'].describe())
    
    # Save output
    adata_tls.write('results/tls_clustered.h5ad')

