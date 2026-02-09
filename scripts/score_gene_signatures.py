"""
Gene Signature Scoring: Calculate and store gene signature scores in AnnData.

This script scores all gene signatures defined in metadata/gene_signatures.json,
including macrophage states (Alveolar macrophages Merad, Mac.2 MoMac M2-like, 
Mac.6 MoMAc M2-like) and other immune/tumor programs.

Uses Scanpy's score_genes function (based on Seurat's AddModuleScore) which:
- Matches control genes by expression level (binned by average expression)
- Calculates score as: mean(signature genes) - mean(matched control genes)
- Controls for technical variation (sequencing depth, batch effects)
- Reduces bias from highly expressed genes dominating the score
- More robust than simple mean averaging

The scores are stored in adata.obs with prefix "sig:" and z-scored versions 
with suffix "_z" for consistent visualization.

Key Advantages over Simple Mean:
1. Expression-level matching: Control genes are sampled from same expression bins
2. Technical variation control: Accounts for sequencing depth and batch effects
3. Reduced bias: Prevents highly expressed genes from dominating scores
4. Standardized approach: Consistent with Seurat/Scanpy best practices
"""

import scanpy as sc
import json
import numpy as np
import pandas as pd
from typing import Mapping, Any, Iterable, Tuple, List, Dict
from pathlib import Path

sc.settings.set_figure_params(dpi=300, dpi_save=400)

# ============================================================================
# Utility Functions for Signature Scoring
# ============================================================================

def _flatten_nested_dict(d: Mapping[str, Any], prefix: Tuple[str, ...] = ()) -> List[Tuple[Tuple[str, ...], List[str]]]:
    """Flatten nested dictionary structure into list of (path, genes) tuples."""
    out = []
    for k, v in d.items():
        if isinstance(v, Mapping):
            out.extend(_flatten_nested_dict(v, prefix + (k,)))
        elif isinstance(v, list):
            out.append((prefix + (k,), [g for g in v if isinstance(g, str)]))
    return out


def _index_genes(genes: Iterable[str], var_names: np.ndarray) -> Tuple[np.ndarray, List[str], List[str]]:
    """Find indices of genes in var_names, return indices, present genes, and missing genes."""
    vm = {g.upper(): i for i, g in enumerate(var_names)}
    idx, present, missing = [], [], []
    for g in genes:
        key = g.upper()
        if key in vm:
            idx.append(vm[key])
            present.append(var_names[vm[key]])
        else:
            missing.append(g)
    return np.array(sorted(set(idx))), present, missing


# Note: _get_matrix and _zscore_cols functions removed as they are no longer needed.
# Scanpy's score_genes handles expression binning and control gene matching internally.


def score_signatures_nested(
    adata,
    signatures_nested: Dict[str, Any],
    prefix: str = "sig",
    use_raw: bool = True,
    ctrl_size: int = 50,
    n_bins: int = 25,
    min_genes: int = 3,
    copy: bool = False,
) -> pd.DataFrame:
    """
    Score nested gene signatures using Scanpy's score_genes (Seurat-based method).
    
    Uses Scanpy's implementation of Seurat's AddModuleScore, which:
    - Bins genes by average expression across cells/spots
    - Samples control genes from same expression bins as signature genes
    - Calculates: score = mean(signature genes) - mean(matched control genes)
    - Controls for technical variation and expression magnitude differences
    
    Parameters
    ----------
    adata : AnnData
        Annotated data object. Must have raw counts if use_raw=True.
    signatures_nested : Dict
        Nested dictionary of gene signatures (e.g., from JSON)
    prefix : str
        Prefix for column names (default: "sig")
    use_raw : bool
        Use adata.raw for expression (default: True).
        score_genes requires raw counts for proper expression binning.
    ctrl_size : int
        Number of control genes to sample per signature gene (default: 50)
    n_bins : int
        Number of expression bins for matching control genes (default: 25)
    min_genes : int
        Minimum number of genes required to score (default: 3)
    copy : bool
        Whether to modify adata in-place or return a copy (default: False)
    
    Returns
    -------
    pd.DataFrame
        DataFrame with signature scores (also stored in adata.obs)
    
    Notes
    -----
    This method is more robust than simple mean averaging because it:
    1. Controls for gene expression magnitude differences
    2. Accounts for technical variation (sequencing depth, batch effects)
    3. Reduces bias from highly expressed genes
    4. Follows Seurat/Scanpy best practices for signature scoring
    
    The algorithm works by:
    - Binning all genes by their average expression across cells/spots
    - For each signature gene, randomly sampling control genes from the same bin
    - Computing per-cell scores as: mean(signature) - mean(controls)
    """
    if copy:
        adata = adata.copy()
    
    leaves = _flatten_nested_dict(signatures_nested)
    assert len(leaves) > 0, "No gene lists in signatures_nested"
    
    # Ensure we have raw counts for proper binning
    if use_raw and adata.raw is None:
        raise ValueError("adata.raw is None but use_raw=True. "
                        "score_genes requires raw counts for expression binning.")
    
    # Get gene names from appropriate source
    if use_raw:
        var_names = np.array(adata.raw.var_names, str)
    else:
        var_names = np.array(adata.var_names, str)
    
    results, report = {}, []
    
    for path, genes in leaves:
        # Find present genes (case-insensitive matching)
        idx, present, missing = _index_genes(genes, var_names)
        col = f"{prefix}:{'/'.join(path)}"
        
        if len(present) < min_genes:
            results[col] = np.full(adata.n_obs, np.nan)
            report.append({
                "column": col, 
                "n_present": len(present), 
                "n_missing": len(missing),
                "status": "skipped"
            })
            continue
        
        # Use Scanpy's score_genes function (Seurat-based)
        # This handles control gene matching and subtraction internally
        try:
            sc.tl.score_genes(
                adata,
                gene_list=present,  # Only use genes present in dataset
                score_name=col,
                ctrl_size=ctrl_size,
                n_bins=n_bins,
                use_raw=use_raw,
            )
            # Extract the score from adata.obs
            if col in adata.obs.columns:
                results[col] = adata.obs[col].values
                report.append({
                    "column": col,
                    "n_present": len(present),
                    "n_missing": len(missing),
                    "status": "ok"
                })
            else:
                # Fallback if score_genes didn't create the column
                results[col] = np.full(adata.n_obs, np.nan)
                report.append({
                    "column": col,
                    "n_present": len(present),
                    "n_missing": len(missing),
                    "status": "error"
                })
        except Exception as e:
            # If score_genes fails (e.g., too few genes, insufficient controls),
            # fall back to NaN
            print(f"Warning: score_genes failed for {col}: {e}")
            results[col] = np.full(adata.n_obs, np.nan)
            report.append({
                "column": col,
                "n_present": len(present),
                "n_missing": len(missing),
                "status": f"error: {str(e)[:50]}"
            })
    
    # Create DataFrame from results (some may already be in adata.obs from score_genes)
    df = pd.DataFrame(results, index=adata.obs_names)
    
    # Ensure all scores are in adata.obs (score_genes already added some)
    for c in df.columns:
        adata.obs[c] = df[c].astype(float)
    
    adata.uns[prefix] = {"report": pd.DataFrame(report)}
    return df


# ============================================================================
# Main Workflow
# ============================================================================

def main():
    """Main function to score gene signatures and save AnnData."""
    print("="*70)
    print("Gene Signature Scoring")
    print("="*70)
    
    # Load gene signatures
    print("\n1. Loading gene signatures...")
    with open('metadata/gene_signatures.json', 'r') as file:
        spatial_signatures = json.load(file)
    print(f"   Loaded {len(spatial_signatures)} signature categories")
    
    # Load AnnData
    print("\n2. Loading AnnData...")
    adata_path = Path('results/adata.annotation.masked.h5ad')
    if not adata_path.exists():
        raise FileNotFoundError(f"AnnData file not found: {adata_path}")
    
    adata = sc.read_h5ad(adata_path)
    adata.var_names_make_unique()
    print(f"   Loaded: {adata.shape} (spots x genes)")
    
    # Score signatures using Scanpy's score_genes (Seurat-based method)
    print("\n3. Scoring gene signatures using Scanpy's score_genes (Seurat-based)...")
    print("   Method: Controls for expression magnitude and technical variation")
    print("   Algorithm: mean(signature genes) - mean(matched control genes)")
    scores = score_signatures_nested(
        adata,
        signatures_nested=spatial_signatures,
        prefix="sig",
        use_raw=True,
        ctrl_size=50,  # Number of control genes per signature gene
        n_bins=25,     # Expression bins for matching
        min_genes=3,
        copy=False,
    )
    print(f"   Scored {len(scores.columns)} signatures")
    
    # Print report
    if "sig" in adata.uns and "report" in adata.uns["sig"]:
        report = adata.uns["sig"]["report"]
        print(f"\n   Signature scoring report:")
        print(f"   - Successfully scored: {len(report[report['status'] == 'ok'])}")
        print(f"   - Skipped (too few genes): {len(report[report['status'] == 'skipped'])}")
        if len(report[report['status'] == 'skipped']) > 0:
            skipped = report[report['status'] == 'skipped'][['column', 'n_present']]
            print(f"\n   Skipped signatures:")
            for _, row in skipped.iterrows():
                print(f"     - {row['column']}: {row['n_present']} genes present")
    
    # Z-score each signature across spots for consistent visualization
    print("\n4. Z-scoring signatures across spots...")
    sig_cols = [c for c in adata.obs.columns if c.startswith("sig:")]
    zscored_count = 0
    for c in sig_cols:
        x = adata.obs[c].values
        if not np.all(np.isnan(x)):
            adata.obs[c + "_z"] = (x - np.nanmean(x)) / (np.nanstd(x) + 1e-8)
            zscored_count += 1
    print(f"   Z-scored {zscored_count} signatures")
    
    # Verify macrophage signatures are present (from Mission.md)
    print("\n5. Verifying macrophage signatures (Mission.md requirements)...")
    required_macrophage_sigs = [
        "sig:Liron/Alveolar macrophages Merad",
        "sig:Liron/Mac.2 MoMac M2-like",
        "sig:Liron/Mac.6 MoMAc M2-like",
    ]
    present_mac = []
    missing_mac = []
    for sig in required_macrophage_sigs:
        if sig in adata.obs.columns:
            present_mac.append(sig)
        else:
            missing_mac.append(sig)
    
    if present_mac:
        print(f"   ✓ Found {len(present_mac)} macrophage signatures:")
        for sig in present_mac:
            print(f"     - {sig}")
    if missing_mac:
        print(f"   ⚠ Missing {len(missing_mac)} macrophage signatures:")
        for sig in missing_mac:
            print(f"     - {sig}")
    
    # Save AnnData with scores
    print("\n6. Saving AnnData with signature scores...")
    output_path = Path('results/adata.img.genescores.h5ad')
    adata.write(output_path)
    print(f"   Saved: {output_path}")
    print(f"   Shape: {adata.shape}")
    print(f"   Signature columns: {len(sig_cols)}")
    print(f"   Z-scored columns: {len([c for c in adata.obs.columns if c.endswith('_z')])}")
    
    print("\n" + "="*70)
    print("✅ Gene signature scoring complete!")
    print("="*70)


if __name__ == "__main__":
    main()

