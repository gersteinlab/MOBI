"""Functions for stats
"""

import numpy as np
from scipy.stats import binom_test

def sig_test(array1, array2):
    """Significant test for two count array1 via bionomial test

    Parameters
    ----------
    array1 : list of int
        array with presumably larger mean
    array2 : list of int
        array with presumably smaller mean

    Returns
    -------
    float
        one-sided p-value
    """
    x = np.sum((array1 - array2) > 0)
    n = np.sum((array1 - array2) != 0)
    pval = binom_test(x, n, alternative='two-sided')

    if (array1.mean() - array2.mean()) > 0:
        pval = pval/2
    else:
        pval = 1
    return pval


def stat_SE_bootstrap(stat):
    """Calculate the standard error of bootstrap estimates of the statistics

    Parameters
    ----------------
    stat : np.array
        list of the statistics from bootstrap, len = bootstrap times

    Returns
    ----------------
    float
        SE
    """
    return np.sqrt(np.sum((stat - stat.mean())**2) / (len(stat)-1))
