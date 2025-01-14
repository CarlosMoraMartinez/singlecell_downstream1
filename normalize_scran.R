library(Seurat)
library(scran)
library(BiocParallel)
library(reticulate)
library(optparse)

scipy_sparse = import("scipy.sparse")

# https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html
option_list <- list(
  make_option(c("-n", "--samplename"), type = "character", default = "",
              help = "Sample name", metavar = "SAMPLENAME"),
  make_option(c("-i", "--input"), type = "character", default = "",
              help = "Input data sparse matrix (.npz)", metavar = "RAWDATA"),
  make_option(c("-o", "--outdir"), type = "character", default = "",
              help = "Output directory", metavar = "OUTDIR"),
  make_option(c("-g", "--groups"), type = "character", default = "",
              help = "Groups data frame", metavar = "GENES")
)
# Create the parser object
opt_parser <- OptionParser(option_list = option_list)
# Parse the arguments
opt <- parse_args(opt_parser)


# Read the input data
samplename <- opt$samplename
infile <- opt$input
groupfile <- opt$groups

outdir <- opt$outdir
if(!dir.exists(outdir)){
  dir.create(outdir, recursive = TRUE)
}

data_mat <- scipy_sparse$load_npz(infile)
input_groups <- read.csv(groupfile, header = TRUE)

size_factors = sizeFactors(
  computeSumFactors(
    SingleCellExperiment(
      list(counts=data_mat)), 
    clusters = input_groups$group,
    min.mean = 0.1,
    BPPARAM = MulticoreParam()
  )
)


resdf = data.frame(size_factors = size_factors)
write.table(resdf, file = paste0(outdir, "/", samplename, "_scran_size_factors.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Scran size factor estimation completed\n")
