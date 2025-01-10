library(tidyverse)
library(SoupX)
library(Seurat)
library(reticulate)
library(optparse)

scipy_sparse = import("scipy.sparse")

# https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html
option_list <- list(
    make_option(c("-n", "--samplename"), type = "character", default = "",
              help = "Sample name", metavar = "SAMPLENAME"),
    make_option(c("-r", "--rawdata"), type = "character", default = "",
              help = "Raw data sparse matrix (.npz)", metavar = "RAWDATA"),
    make_option(c("-p", "--preprocessed"), type = "character", default = "",
              help = "Preprocessed data fromm preprocess_1 script, sparse matrix (.npz)", metavar = "PREPDATA"),
    make_option(c("-c", "--cells"), type = "character", default = "",
              help = "Cells data frame", metavar = "CELLS"),
    make_option(c("-g", "--genes"), type = "character", default = "",
              help = "Genes data frame", metavar = "GENES"),
    make_option(c("-x", "--clusters"), type = "character", default = "",
              help = "Cluster per cell", metavar = "CLUSTERS"),
    make_option(c("-o", "--outdir"), type = "character", default = "",
              help = "Output directory", metavar = "OUTDIR")
)
# Create the parser object
opt_parser <- OptionParser(option_list = option_list)
# Parse the arguments
opt <- parse_args(opt_parser)


# Read the input data
samplename <- opt$samplename
outdir <- opt$outdir
if(!dir.exists(outdir)){
    dir.create(outdir, recursive = TRUE)
}

data_tod <- scipy_sparse$load_npz(opt$rawdata)
data <- scipy_sparse$load_npz(opt$preprocessed)
cells <- read.csv(opt$cells, header = TRUE)
genes <- read.csv(opt$genes, header = TRUE)
soupx_groups <- read.csv(opt$clusters, header = TRUE)

# specify row and column names of data
rownames(data) = genes$gene
colnames(data) = cells$cell
# ensure correct sparse format for table of counts and table of droplets
data <- as(data, "sparseMatrix")
data_tod <- as(data_tod, "sparseMatrix")

# Generate SoupChannel Object for SoupX 
sc = SoupChannel(data_tod, data, calcSoupProfile = FALSE)

# Add extra meta data to the SoupChannel object
soupProf = data.frame(row.names = rownames(data), est = rowSums(data)/sum(data), counts = rowSums(data))
sc = setSoupProfile(sc, soupProf)
# Set cluster information in SoupChannel
sc = setClusters(sc, soupx_groups$soupx_groups)

# Estimate contamination fraction
sc  = autoEstCont(sc, doPlot=FALSE)
# Infer corrected table of counts and rount to integer
out = adjustCounts(sc, roundToInt = TRUE)

scipy_sparse$save_npz(paste0(outdir, "/", samplename, "_soupX.npz"), out)

cat(paste0("Output saved to ", paste0(outdir, "/", samplename, "_soupX.npz"), "\n"))