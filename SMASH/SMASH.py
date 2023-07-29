''' Function for the main test
'''

import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from Kernel_definitions import get_l_limits, All_kernel_R2
from Pvalue_adjustments import ACAT, Min_p
from Testing_functions import SMASH_linear, SMASH_linear_Gauss, SMASH_linear_Cosine, SMASH_single_asym_test

def SMASH(Y, Cords, len_l = 10, mean_only = "False", kernel_covariance = "All", forcePD = "False"):
     ''' Run SMASH on the gene and location dataframes
    
        Parameters
        ----------
        Y: A pandas dataframe of the gene expression profiles, rows representing cells/spots and columns represeneting genes
        
        Cords: A pandas dataframe of 2D spatial co-ordinates (xy co-ordinates), rownames of Cords should match the rownames of Y
                       
        len_l: An integer specifying how many lengthscale and period values to consider in the kernel covariance matrices.
               Default is 10
        
        mean_only: A binary string that shall be set to be "FALSE" to take into account Gaussian and Cosine kernel 
                   covariance matrices. If set to be "TRUE", only linear kernel i.e., mean-based test will be performed
        
        kernel_covariance: A string specifying if both Gaussian and Cosine kernel covariance matrices to be used. Default is "All"
                which implies both the types of kernel covariance matrices will be used. If mean_only = "TRUE", this
                parameter is redundant as neither of the kernel covariance matrices will be used
        
        forcePD: "True" or" "False" where the latter implies that the cosine kernel matrix will be considered as is, 
                 without any adjustement to force it to be positive definite. If "True", eigenvalue-truncation based adjustment                    
                 will be performed.        
                    
        Returns
        -------
        results: A table of spatial statistics for each gene.
     '''
     Gene_names = Y.columns
     Y = Y.values
     Cords = Cords.values
     
     # Get upper and lower limits of the lengthscale parameter and the distance matrix 
     lp_min, lp_max, R2 = get_l_limits(Cords)
     
     # Choose ten (or length_l many) uniformly varying values of the lengthscale between the limits
     l_series =  np.logspace(np.log10(lp_min), np.log10(lp_max), len_l) 
     Cords = Cords.astype("float")
     
     # Mean and SD scaling of the gene expression data
     scaler = StandardScaler(with_mean = True, with_std = True)
     Y_scaled = scaler.fit_transform(Y)
     p = Y_scaled.shape[1]; n = Y_scaled.shape[0]
     
     # Perform the mean-based (first order) test (eqv. to SPARK-X)
     print("Starting Mean Tests!")
     lin_direct_p = SMASH_linear(Y_scaled, Cords) # Test using linear kernel directly on the co-ordinates
     lin_gauss_transformed_p = SMASH_linear_Gauss(Y_scaled, Cords) # Test using linear kernel on five Gaussian transformed versions of the co-ordinates
     lin_cosine_transformed_p = SMASH_linear_Cosine(Y_scaled, Cords) # Test using linear kernel on five Cosine transformed versions of the co-ordinates
     lin_gauss_cosine_p = ACAT(np.vstack((lin_gauss_transformed_p, lin_cosine_transformed_p))) # Combine the last two tests (their p-values) using a Cauchy combination rule
     lin_all_combined_p = ACAT(np.vstack((lin_direct_p, lin_gauss_transformed_p, lin_cosine_transformed_p))) # Combine all three tests (their p-values) using a Cauchy combination rule (eqv. to SPARK-X)
     print("Finished Mean Tests!")
     
     # Perform the covariance-based test involving Gaussian and Cosine kernel covariance matrices 
     if mean_only == "True": # "TRUE" if only a mean-based test (SPARK-X) is sought, not recommended due to having low power 
         SPARKX = pd.DataFrame(np.hstack((np.array(Gene_names).reshape(p, 1), lin_all_combined_p.reshape(p, 1)))); SPARKX = SPARKX.set_axis(['Gene', 'p-val'], axis = 1)
         return {'SPARK-X': lin_all_combined_p}
     else:
         print("Starting Covariance Tests!")
         if kernel_covariance == "All": # "All" implies that both Gaussian and Cosine kernel covariance matrices will be used
            pvals = np.zeros(( 2*len(l_series), p)); r = 0; 
            for l in tqdm(l_series): # Loop over ten (or, len_l many) values of the lengthscale (period)
              dist_both =  All_kernel_R2(R2, l, which = "Both", forcePD = forcePD) 
              dist1 = dist_both[0]; dist2 = dist_both[1] # Obtain the Gaussian (at index 0) and Cosine (at index 1) kernel covariance matrices for a particular value of the lengthscale (or, period) l
              dist1 = dist1 - 1/n*np.repeat(np.sum(dist1, axis = 0)[np.newaxis,:], n, axis = 0) # Row and column scaling of the Gaussian kernel covariance matrix
              dist1 = dist1 - 1/n*np.repeat(np.sum(dist1, axis = 1)[:,np.newaxis], n, axis = 1) 
              dist2 = dist2 - 1/n*np.repeat(np.sum(dist2, axis = 0)[np.newaxis,:], n, axis = 0) # Row and column scaling of the Gaussian kernel covariance matrix   
              dist2 = dist2 - 1/n*np.repeat(np.sum(dist2, axis = 1)[:,np.newaxis], n, axis = 1)          
              trace_dist = np.sum(np.diag(dist1)); trace_dist2 = np.sum(dist1 *  dist1) # Obtain the parameters of the approximate asymptotic Gamma distribution of the test statistic with the Gaussuan kernel under null
              k1 = trace_dist**2/2/trace_dist2; theta1 = 2*trace_dist2/trace_dist
              trace_dist = np.sum(np.diag(dist2)); trace_dist2 = np.sum(dist2 *  dist2) # Obtain the parameters of the approximate asymptotic Gamma distribution of the test statistic with the Cosine kernel under null      
              k2 = trace_dist**2/2/trace_dist2; theta2 = 2*trace_dist2/trace_dist
              pvals[2*r,:] = SMASH_single_asym_test(Y_scaled, dist1, k1, theta1); # Perform the test and obtain the p-values for the Gaussian kernel
              pvals[2*r+1,:] = SMASH_single_asym_test(Y_scaled, dist2, k2, theta2); # Perform the test and obtain the p-values for the Cosine kernel
              r = r + 1
         elif kernel_covariance == "Gaussian": # "Gaussian" implies that only Gaussian kernel covariance matrices will be used
          pvals = np.zeros(( len(l_series), p)); r = 0; 
          for l in tqdm(l_series):
              dist1 =  All_kernel_R2(R2, l, which = "Gaussian")[0]
              dist1 = dist1 - 1/n*np.repeat(np.sum(dist1, axis = 0)[np.newaxis,:], n, axis = 0)
              dist1 = dist1 - 1/n*np.repeat(np.sum(dist1, axis = 1)[:,np.newaxis], n, axis = 1)
              trace_dist = np.sum(np.diag(dist1)); trace_dist2 = np.sum(dist1 *  dist1)
              k1 = trace_dist**2/2/trace_dist2; theta1 = 2*trace_dist2/trace_dist
              pvals[r,:] = SMASH_single_asym_test(Y_scaled, dist1, k1, theta1);
              r = r + 1
         elif kernel_covariance == "Cosine":  # "Cosine" implies that only Cosine kernel covariance matrices will be used
          pvals = np.zeros(( len(l_series), p)); r = 0; 
          for l in tqdm(l_series):
              dist1 =  All_kernel_R2(R2, l, which = "Cosine", forcePD = forcePD)[0]
              dist1 = dist1 - 1/n*np.repeat(np.sum(dist1, axis = 0)[np.newaxis,:], n, axis = 0)
              dist1 = dist1 - 1/n*np.repeat(np.sum(dist1, axis = 1)[:,np.newaxis], n, axis = 1)
              trace_dist = np.sum(np.diag(dist1)); trace_dist2 = np.sum(dist1 *  dist1)
              k1 = trace_dist**2/2/trace_dist2; theta1 = 2*trace_dist2/trace_dist
              pvals[r,:] = SMASH_single_asym_test(Y_scaled, dist1, k1, theta1);
              r = r + 1
         comb_pvals = Min_p(np.vstack((lin_direct_p.reshape(1, p), pvals))) # Combine three results (p-values) using minimum p-value approach from 
                                                                            # 1) the mean-based test using linear kernel directly on co-ordinates
                                                                            # 2) the covariance-based test using Gaussian kernel 
                                                                            # 3) the covariance-based test using Cosine kernel
                                                                            # This p-value can be roughly interpreted as the one we get from SpatialDE
         final_pvals = Min_p(np.vstack((lin_gauss_cosine_p, comb_pvals)))   # Finally combine the above p-values with the two other linear kernel based p-values
         print("Finished Covariance Tests!")    
         SPARKX = pd.DataFrame(np.hstack((np.array(Gene_names).reshape(p, 1), lin_all_combined_p.reshape(p, 1)))); SPARKX = SPARKX.set_axis(['Gene', 'p-val'], axis = 1)
         SpatialDE_approx = pd.DataFrame(np.hstack((np.array(Gene_names).reshape(p, 1), comb_pvals.reshape(p, 1)))); SpatialDE_approx = SpatialDE_approx.set_axis(['Gene', 'p-val'], axis = 1)
         SMASH = pd.DataFrame(np.hstack((np.array(Gene_names).reshape(p, 1), final_pvals.reshape(p, 1)))); SMASH = SMASH.set_axis(['Gene', 'p-val'], axis = 1)
         return {'SPARK-X Approx.': SPARKX, 'SpatialDE Approx.': SpatialDE_approx, 'SMASH': SMASH}
        

def Expression_plot(Data, Gene_name, s = 0.5, cmap = 'viridis_r'):
    plt.figure()
    ax = plt.gca()
    ax.invert_yaxis()
    ax.set_title(label = Gene_name, loc = 'center')
    im = ax.scatter(Data['X'], Data['Y'], s = s, alpha = 1, c = Data[Gene_name]/max(Data[Gene_name]), cmap = cmap)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)     
    plt.colorbar(im, cax=cax)
    return None
