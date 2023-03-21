# Python implementation of SMASH

### How to install?
We recommend downloading the Github repository as a ZIP file and unpacking it. The users need to change the path (which currently looks like, path = "/Users/sealso/Documents/GitHub/SMASH-package") to the location of the unpacked SMASH-package folder in their system. 


### Overview of the main functions
-  The function SMASH fits our proposed method. Along with the p-values corresponding to the method SMASH, the function also returns two additional results, one corresponding to the method SPARK-X and the other result corresponding to an approximate version of the method SpatialDE. 

-  The function Expression_plot can be used for basic visualization of a gene expression  in cells/spots.


### Jupyter notebooks with function usage and two real data analyses
- The Jupyter notebook entitled "Merfish_analysis.ipynb" provides a thorough guide on how to use the package on a mouse cerebellum data collected using the MERFISH platform. This notebook explains all the arguments that can be tweaked in the main functions. 

- The Jupyter notebook entitled "Visium_analysis.ipynb" provides a guide on how to use the package on a human DLPFC data collected using the 10X Visium platform. 


### Python modules required to be installed
- The package requires the following module to be pre-installed,
  1. matplotlib, install using: "conda install -c conda-forge matplotlib"  (https://matplotlib.org/)

* We recommend using Anaconda (https://www.anaconda.com/products/individual) and Python version > 3.9. 


### References

1. Zhang, K., Peters, J., Janzing, D., Sch¨olkopf, B.: Kernel-based conditional independence test and application in causal discovery. arXiv preprint arXiv:1202.3775 (2012)

2. Maynard, K.R., Collado-Torres, L., Weber, L.M., Uytingco, C., Barry, B.K., Williams, S.R., Catallini, J.L., Tran, M.N., Besich, Z., Tippani, M., et al.: Transcriptome-scale spatial gene expression in the human dorsolateral
prefrontal cortex. Nature neuroscience 24(3), 425–436 (2021)

3. Moffitt, J.R., Bambah-Mukku, D., Eichhorn, S.W., Vaughn, E., Shekhar, K., Perez, J.D., Rubinstein, N.D., Hao, J., Regev, A., Dulac, C., et al.: Molecular, spatial, and functional single-cell profiling of the hypothalamic
preoptic region. Science 362(6416), 5324 (2018)

4. Zhu, J., Sun, S., Zhou, X.: Spark-x: non-parametric modeling enables scalable and robust detection of spatial expression patterns for large spatial transcriptomic studies. Genome Biology 22(1), 1–25 (2021)

5. Svensson, V., Teichmann, S.A., Stegle, O.: Spatialde: identification of spatially variable genes. Nature methods 15(5), 343–346 (2018)

