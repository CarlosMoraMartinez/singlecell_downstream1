library(tidyverse)
library(Seurat)

library(SeuratDisk)

outdir <- "/home/carlos/projects/singlecell_Dec2024/secondary_analysis/cellranger1_secondary1/"
samples <- c("A", "B", "C")

sample <- "A"

infile <- paste0(outdir, sample, "/", sample, "_soupX_MarkedDoublets_Norm.h5ad")
outfile <- paste0(outdir, sample, "/", sample, "_soupX_MarkedDoublets_Norm.h5seurat")
Convert(infile, outfile, assay="RNA", overwrite = TRUE)
# This creates a copy of this .h5ad object reformatted into .h5seurat inside the example_dir directory

# This .d5seurat object can then be read in manually
seuratObject <- LoadH5Seurat(paste0(outdir, sample, "/", sample, "_soupX_MarkedDoublets_Norm.h5seurat"))
