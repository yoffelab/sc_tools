"""
Statistical testing utilities.

Provides functions for:
- Mann-Whitney U test (mwu)
- FDR correction (Benjamini-Hochberg)
- Other statistical tests
"""

from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from typing import List, Tuple, Optional
import numpy as np
import pandas as pd


def mwu(data1: np.ndarray, data2: np.ndarray, 
        alternative: str = 'two-sided') -> Tuple[float, float]:
    """
    Perform Mann-Whitney U test.
    
    Parameters
    ----------
    data1 : array-like
        First sample
    data2 : array-like
        Second sample
    alternative : str
        Alternative hypothesis: 'two-sided', 'less', or 'greater'
    
    Returns
    -------
    tuple
        (statistic, pvalue)
    """
    statistic, pvalue = mannwhitneyu(data1, data2, alternative=alternative)
    return statistic, pvalue


def fdr_correction(pvalues: np.ndarray, 
                   method: str = 'fdr_bh',
                   alpha: float = 0.05) -> Tuple[np.ndarray, np.ndarray, float, float]:
    """
    Apply FDR (False Discovery Rate) correction using Benjamini-Hochberg method.
    
    Parameters
    ----------
    pvalues : array-like
        Array of p-values to correct
    method : str
        Correction method ('fdr_bh' for Benjamini-Hochberg)
    alpha : float
        Significance level
    
    Returns
    -------
    tuple
        (rejected, pvals_corrected, alpha_sidak, alpha_bonf)
    """
    rejected, pvals_corrected, alpha_sidak, alpha_bonf = multipletests(
        pvalues, alpha=alpha, method=method
    )
    return rejected, pvals_corrected, alpha_sidak, alpha_bonf


__all__ = ['mwu', 'fdr_correction']
