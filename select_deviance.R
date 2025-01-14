library(Seurat)
library(scran)
library(BiocParallel)
library(reticulate)
library(optparse)
library(scry) # feature selection

scipy_sparse = import("scipy.sparse")

# https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html
option_list <- list(
  make_option(c("-n", "--samplename"), type = "character", default = "",
              help = "Sample name", metavar = "SAMPLENAME"),
  make_option(c("-i", "--input"), type = "character", default = "",
              help = "Input data sparse matrix (.npz)", metavar = "RAWDATA"),
  make_option(c("-o", "--outdir"), type = "character", default = "",
              help = "Output directory", metavar = "OUTDIR")
)
# Create the parser object
opt_parser <- OptionParser(option_list = option_list)
# Parse the arguments
opt <- parse_args(opt_parser)


# Read the input data
samplename <- opt$samplename
infile <- opt$input

outdir <- opt$outdir
if(!dir.exists(outdir)){
  dir.create(outdir, recursive = TRUE)
}

data_mat <- scipy_sparse$load_npz(infile)


## Find genes with greater deviance
sce = devianceFeatureSelection(SingleCellExperiment(
  list(counts=data_mat)))

write.table(rowData(sce), file = paste0(outdir, "/", samplename, "_deviance_genes.tsv"), 
            sep = "\t", quote = FALSE, row.names = TRUE)
