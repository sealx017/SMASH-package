''' Functions to combine multiple p-values
'''

import numpy as np
from scipy.stats import cauchy


def ACAT(Pvals):
    '''
    Returns combined p-value following Cauchy combination rule, used only for linear kernel-based tests
    
    '''
    Weights = 1/Pvals.shape[0]; p = Pvals.shape[1];
    Cauchy_stat = np.zeros(p);
    Cauchy_pval = np.zeros(p);
    Pvals[Pvals == 0] = np.min(Pvals[Pvals != 0])/2
    Pvals[Pvals < 1e-250] = 1e-250
    for r in range(p):
     is_small = Pvals[:,r] < 1e-16; 
     if (np.sum(is_small)==0):
        Cauchy_stat[r]  = np.sum(Weights*np.tan((0.5-Pvals[:,r])*np.pi))
     else:
        Cauchy_stat[r] = np.sum((Weights/Pvals[is_small,r])/np.pi) + np.sum(Weights*np.tan((0.5-Pvals[~ is_small,r])*np.pi))
     if Cauchy_stat[r] > 1e+14:
        Cauchy_pval[r] = 1/Cauchy_stat[r]/np.pi
     else:
        Cauchy_pval[r] = np.exp(cauchy.logsf(Cauchy_stat[r])) 
    return(Cauchy_pval)


def Min_p(Pvals):
    '''
    Returns combined p-value following  minimum p-value rule, used for non-linear kernel-based tests and the final p-value
    
    '''
    p = Pvals.shape[1]; k = Pvals.shape[0]
    min_pval = np.zeros(p);
    for r in range(p):
     min_pval[r] = min(np.min(k*np.min(Pvals[:,r])), 1)
    return(min_pval)