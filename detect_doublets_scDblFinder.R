library(Seurat)
library(scater)
library(scDblFinder)
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
  make_option(c("-s", "--seed"), type = "character", default = "",
              help = "Random seed", metavar = "SEED")
)
# Create the parser object
opt_parser <- OptionParser(option_list = option_list)
# Parse the arguments
opt <- parse_args(opt_parser)


# Read the input data
samplename <- opt$samplename
infile <- opt$input
seed <- opt$seed

outdir <- opt$outdir
if(!dir.exists(outdir)){
  dir.create(outdir, recursive = TRUE)
}

data_mat <- scipy_sparse$load_npz(infile)

set.seed(seed)
sce = scDblFinder(
  SingleCellExperiment(
    list(counts=data_mat),
  ) 
)
doublet_score = sce$scDblFinder.score
doublet_class = sce$scDblFinder.class

tab <- data.frame(table(doublet_class) )
write.table(tab, file = paste0(outdir, "/", samplename, "_doublet_class_table.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

resdf = data.frame(doublet_score = doublet_score, doublet_class = doublet_class)
write.table(resdf, file = paste0(outdir, "/", samplename, "_doublet_score.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Doublet detection completed\n")
