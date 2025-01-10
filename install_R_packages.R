install.packages("devtools")
install.packages("tidyverse")
install.packages("BiocManager")
install.packages("hdf5r")
install.packages('Seurat')
install.packages("sctransform")
install.packages("remotes")
install.packages("optparse")
remotes::install_github("10xGenomics/loupeR")
loupeR::setup()


install.packages('SoupX')
BiocManager::install("scDblFinder")
BiocManager::install('glmGamPoi')
BiocManager::install('scran')
remotes::install_github("mojaveazure/seurat-disk")
devtools::install_github('satijalab/seurat-data')

BiocManager::install('zellkonverter')
# sce <- zellkonverter::readH5AD(h5ad_file, verbose = TRUE) # from scanpy HD5F file

install.packages('reticulate')
library(reticulate)

py_install("pandas")
py_install("scipy")