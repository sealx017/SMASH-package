# Python implementation of SMASH

### How to install?
We recommend downloading the Github repository as a ZIP file and unpacking it. The users need to change the "path" (currently, path = "/Users/sealso/Documents/GitHub/SMASH-package") mentioned at the start of each notebook  to the location of the unpacked SMASH-package-main folder in their system. 


### Overview of the main functions
-  The function SMASH fits our proposed method. Along with the p-values corresponding to SMASH, the function also returns two additional results corresponding to approximate versions of the two methods SPARK-X [1] and SpatialDE [2]. 

-  The function Expression_plot can be used for basic visualization of a gene expression in cells/spots.


### Jupyter notebooks with function usage and two real data analyses
- The Jupyter notebook entitled [Merfish_analysis.ipynb](Merfish_analysis.ipynb) provides a thorough guide on how to use the package on a mouse cerebellum data collected using the MERFISH platform [3]. This notebook explains all the arguments that can be tweaked in the main functions. 

- The Jupyter notebook entitled [Human_DLPFC_Visium_analysis.ipynb](Human_DLPFC_Visium_analysis.ipynb) provides a guide on how to use the package on a human DLPFC data collected using the 10X Visium platform [4]. 


### Python modules required to be installed
- The package and the notebooks provided, requires the following modules to be pre-installed,
  1. matplotlib, install using: "pip install matplotlib"  (https://matplotlib.org/)
  2. matplotlib_venn, install using: "pip install matplotlib-venn"  (https://pypi.org/project/matplotlib-venn/)
  3. blosc, install using: "pip install blosc" (https://anaconda.org/anaconda/blosc)

* We recommend using Anaconda (https://www.anaconda.com/products/individual) and Python version > 3.9. 

### Data provided

Two datasets: (1) a mouse cerebellum data collected using the MERFISH platform [3], and (2)  a human DLPFC data collected using the 10X Visium platform [4], are provided as pickled python objects in the [Data](Data) folder. 

### References

1. Zhu, J., Sun, S., Zhou, X.: Spark-x: non-parametric modeling enables scalable and robust detection of spatial expression patterns for large spatial transcriptomic studies. Genome Biology 22(1), 1–25 (2021)

2. Svensson, V., Teichmann, S.A., Stegle, O.: Spatialde: identification of spatially variable genes. Nature methods 15(5), 343–346 (2018)

3. Maynard, K.R., Collado-Torres, L., Weber, L.M., Uytingco, C., Barry, B.K., Williams, S.R., Catallini, J.L., Tran, M.N., Besich, Z., Tippani, M., et al.: Transcriptome-scale spatial gene expression in the human dorsolateral
prefrontal cortex. Nature neuroscience 24(3), 425–436 (2021)

4. Moffitt, J.R., Bambah-Mukku, D., Eichhorn, S.W., Vaughn, E., Shekhar, K., Perez, J.D., Rubinstein, N.D., Hao, J., Regev, A., Dulac, C., et al.: Molecular, spatial, and functional single-cell profiling of the hypothalamic
preoptic region. Science 362(6416), 5324 (2018)


