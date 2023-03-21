''' Functions for the association tests
'''

import numpy as np
from math import pi
from scipy.stats import chi2, gamma
from scipy.linalg.blas import sgemm
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler
from .Kernel_definitions import get_l_limits, All_kernel_R2
from .Pvalue_adjustments import ACAT, Min_p


def SMASH_linear(Y, Cords):
    p = Y.shape[1]; test = np.zeros(p)
    Cords = Cords - np.mean(Cords, axis = 0)
    X = Cords
    XTX = sgemm(alpha=1, a = X, b =  X, trans_a = 1)
    XTX_inv = np.linalg.inv(XTX)   
    beta = np.dot(XTX_inv, sgemm(alpha=1, a = X, b = Y, trans_a = 1))
    var_beta = XTX
    var_betaTbeta = sgemm(alpha=1,  a = var_beta, b = beta, trans_a = 1)
    for i in range(p):
        test[i] = sgemm(alpha=1, a = beta[:,i], b = var_betaTbeta[:,i], trans_a =1)
    pval = np.exp(chi2.logsf(test, 2))
    return pval.reshape(1, p)


def SMASH_linear_Gauss(Y, Cords):
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

def SMASH_single(Y1, dist):
    p = Y1.shape[1]; statistic = np.zeros(p);
    yTsig= sgemm(alpha=1, a = Y1,
                b =  dist, trans_a = 1)
    for i in range(p):
     statistic[i] = (yTsig[i,:] @ Y1[:,i])
    return statistic
def SPARKY(Y, dist):
    p = Y.shape[1]; statistic = np.zeros(p); 
    yTsig= sgemm(alpha=1, a = Y, b =  dist, trans_a = 1)
    for i in range(p):
     statistic[i] = (yTsig[i,:] @ Y[:,i])  
    return statistic

def SPARKY_single_asym_test(Y, dist, k, theta):
    pval = np.exp(gamma.logsf(SPARKY(Y, dist), k, 0, theta))
    return(pval)


def SMASH_single_asym_test(Y, dist, a, b):
    scaler = StandardScaler()
    Y = scaler.fit_transform(Y)
    pval =  1 - chi2.cdf(SMASH_single(Y, dist)/a, b)
    return(pval)

def SPARKY_single_diff_kernels(Y, Cords, len_l = 10, method = "Gaussian", transform = "True", mean_only = "False", subset = 2500):
     
    ''' Run the SpatialDE test on an AnnData object
    
        Parameters
        ----------
        adata: An AnnData object with counts in the .X field.
        
        coord_columns: A list with the columns of adata.obs which represent spatial
                       coordinates. Default ['x', 'y']
                       
        regress_formula: A patsy formula for linearly regressing out fixed effects
                         from columns in adata.obs before fitting the SpatialDE models.
                         Default is 'np.log(total_counts)'
        Returns
        -------
        results: A table of spatial statistics for each gene.
   '''
     lp_min, lp_max, R2 = get_l_limits(Cords)
     l_series =  np.logspace(np.log10(lp_min), np.log10(lp_max), 10) 
     Cords = Cords.astype("float")
     scaler = StandardScaler(with_mean = True, with_std= True)
     Y_scaled = scaler.fit_transform(Y)
     if subset>0:
         Y_scaled = Y_scaled[:,range(subset)]
     p = Y_scaled.shape[1]; n = Y_scaled.shape[0]
     print("Starting Mean Tests!")
     lin_direct_p = SMASH_linear(Y_scaled, Cords)
     lin_gauss_transformed_p = SMASH_linear_Gauss(Y_scaled, Cords)
     lin_cosine_transformed_p = SMASH_linear_Cosine(Y_scaled, Cords)
     lin_p = ACAT(np.vstack((lin_gauss_transformed_p, lin_cosine_transformed_p))) 
     print("Finished Mean Tests!")
     print("Starting Covariance Tests!")
     if mean_only == "True":
         final_pvals = lin_p
     else:    
         if method == "All":
            pvals = np.zeros(( 2*len(l_series), p)); r = 0; 
            for l in tqdm(l_series):
              dist_both =  All_kernel_R2(R2, l, which = "Both")
              dist1 = dist_both[0]; dist2 = dist_both[1]
              dist1 = dist1 - 1/n*np.repeat(np.sum(dist1, axis = 0)[np.newaxis,:], n, axis = 0)
              dist1 = dist1 - 1/n*np.repeat(np.sum(dist1, axis = 1)[:,np.newaxis], n, axis = 1)
              dist2 = dist2 - 1/n*np.repeat(np.sum(dist2, axis = 0)[np.newaxis,:], n, axis = 0)
              dist2 = dist2 - 1/n*np.repeat(np.sum(dist2, axis = 1)[:,np.newaxis], n, axis = 1)             
              trace_dist = np.sum(np.diag(dist1)); trace_dist2 = np.sum(dist1 *  dist1)
              k1 = trace_dist**2/2/trace_dist2; theta1 = 2*trace_dist2/trace_dist
              trace_dist = np.sum(np.diag(dist2)); trace_dist2 = np.sum(dist2 *  dist2)      
              k2 = trace_dist**2/2/trace_dist2; theta2 = 2*trace_dist2/trace_dist
              pvals[2*r,:] = SPARKY_single_asym_test(Y_scaled, dist1, k1, theta1);
              pvals[2*r+1,:] = SPARKY_single_asym_test(Y_scaled, dist2, k2, theta2); 
              r = r + 1
         elif method == "Gaussian":
          pvals = np.zeros(( len(l_series), p)); r = 0; 
          for l in tqdm(l_series):
              dist1 =  All_kernel_R2(R2, l, which = "Gaussian")[0]
              dist1 = dist1 - 1/n*np.repeat(np.sum(dist1, axis = 0)[np.newaxis,:], n, axis = 0)
              dist1 = dist1 - 1/n*np.repeat(np.sum(dist1, axis = 1)[:,np.newaxis], n, axis = 1)
              trace_dist = np.sum(np.diag(dist1)); trace_dist2 = np.sum(dist1 *  dist1)
              k1 = trace_dist**2/2/trace_dist2; theta1 = 2*trace_dist2/trace_dist
              pvals[r,:] = SPARKY_single_asym_test(Y_scaled, dist1, k1, theta1);
              r = r + 1
         comb_pvals = Min_p(np.vstack((lin_direct_p.reshape(1, p), pvals)))
         final_pvals = np.zeros(p)
         final_pvals = Min_p(np.vstack((lin_p, comb_pvals)))
     print("Finished Covariance Tests!")    
     return [lin_p, comb_pvals, final_pvals]
        


