''' Functions for the individual association tests
'''

import numpy as np
from math import pi
from scipy.stats import chi2, gamma
from scipy.linalg.blas import sgemm


def SMASH_linear(Y, Cords):
    '''
    Returns an array of p-values corresponding to linear kernel based association test between Y and Cords
    
    '''
    p = Y.shape[1]; test = np.zeros(p)
    Cords = Cords - np.mean(Cords, axis = 0)
    X = Cords
    XTX = sgemm(alpha=1, a = X, b =  X, trans_a = 1)
    XTX_inv = np.linalg.inv(XTX)   
    beta = np.dot(XTX_inv, sgemm(alpha=1, a = X, b = Y, trans_a = 1))
    var_beta = XTX/np.mean((Y-np.dot(X, beta))**2)
    var_betaTbeta = sgemm(alpha=1,  a = var_beta, b = beta, trans_a = 1)
    for i in range(p):
        test[i] = sgemm(alpha=1, a = beta[:,i], b = var_betaTbeta[:,i], trans_a =1)
    pval = np.exp(chi2.logsf(test, 2))
    return pval.reshape(1, p)


def SMASH_linear_Gauss(Y, Cords):
    '''
    Returns an array of p-values corresponding to linear kernel based association test between Y and g1(Cords/k)
    where g1 denotes a co-ordinate wise Gaussian transformation and k is a scaling constant. Five data-derived values of k are
    considered
    
    '''    
    p = Y.shape[1];  
    quantile_sequence = np.linspace(0.2,1,5)
    Cords = Cords - np.mean(Cords, axis = 0)
    limits_x = np.quantile(abs(Cords[:,0]), quantile_sequence)
    limits_y = np.quantile(abs(Cords[:,1]), quantile_sequence)
    pvals = np.zeros((5, p))
    for i in range(5):
        Trans_Cords = np.vstack((np.exp(-Cords[:,0]**2/2/limits_x[i]**2), 
        np.exp(-Cords[:,1]**2/2/limits_y[i]**2))).T
        pvals[i,:] = SMASH_linear(Y, Trans_Cords)
    return pvals

def SMASH_linear_Cosine(Y, Cords):
    '''
    Returns an array of p-values corresponding to linear kernel based association test between Y and g2(Cords/k)
    where g2 denotes a co-ordinate wise Cosine transformation and k is a scaling constant. Five data-derived values of k are
    considered
    
    ''' 
    p = Y.shape[1];  
    quantile_sequence = np.linspace(0.2,1,5)
    Cords = Cords - np.mean(Cords, axis = 0)
    limits_x = np.quantile(abs(Cords[:,0]), quantile_sequence)
    limits_y = np.quantile(abs(Cords[:,1]), quantile_sequence)
    pvals = np.zeros((5, p))
    for i in range(5):
        Trans_Cords = np.vstack((np.cos(2*pi*Cords[:,0]/limits_x[i]), 
        np.cos(2*pi*Cords[:,1]/limits_y[i]))).T
        pvals[i,:] = SMASH_linear(Y, Trans_Cords)
    return pvals


def SMASH_single_asym_test(Y, dist, k, theta):
    '''
    Parameters
    ----------
        Y: A numpy array of the gene expression profiles provided from the main function
        dist: A kernel covariance matrix between spatial co-ordinates provided from the main function
        k: The first distributional parameter, precomputed from dist and provided from the main function
        theta: The second distributional parameter, precomputed from dist and provided from the main function
    
    Returns
    -------
        pval: An array of p-values corresponding to an association test between Y and dist
    
    ''' 
    p = Y.shape[1]; statistic = np.zeros(p); 
    yTsig= sgemm(alpha=1, a = Y, b =  dist, trans_a = 1)
    for i in range(p):
     statistic[i] = (yTsig[i,:] @ Y[:,i])  
    pval = np.exp(gamma.logsf(statistic, k, 0, theta))
    return pval

