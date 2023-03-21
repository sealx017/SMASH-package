''' Functions to compute kernel covariance matrices from the spatial co-ordinates
'''

import numpy as np


def get_l_limits(X):
    '''
    Returns the two boundary values of the lengthscale parameter: l_min and l_max, and the distance matrix R2
    
    '''
    Xsq = np.sum(np.square(X), 1)
    R2 = -2. * np.dot(X, X.T) + (Xsq[:, None] + Xsq[None, :])
    R2 = np.clip(R2, 0, np.inf)
    R_vals = np.unique(R2.flatten())
    R_vals = R_vals[R_vals > 1e-8]
    l_min = np.sqrt(R_vals.min()) / 2
    l_max = np.sqrt(R_vals.max()) * 2 
    return l_min, l_max, R2


def All_kernel_R2(R2, l, which = "Both"):
    '''
    Returns Gaussian and Cosine kernel covariance matrices for a particular choice of the lengthscale (or, period) paramter l. 
    which = "Both" returns both the matrices, while the other two options ("Gaussian" or "Cosine") return only the corresponding matrices
    
    '''
    Both_kernels = []
    R2 = R2.astype(float)
    if which == "Both":
      Both_kernels.append(np.exp(-R2 / (2 * l ** 2)))
      Both_kernels.append(np.cos(2 * np.pi * np.sqrt(R2) / l))
    elif which == "Gaussian":
      Both_kernels.append(np.exp(-R2 / (2 * l ** 2)))
    elif which == "Cosine":
      Both_kernels.append(np.cos(2 * np.pi * np.sqrt(R2) / l))
    return Both_kernels