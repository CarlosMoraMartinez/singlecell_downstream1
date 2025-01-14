import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import os
from myutils import *
import seaborn as sns
import scipy.sparse
from scipy.sparse import csr_matrix, issparse
import pandas as pd


MIT_STR = "mt-"
RIBOSOMAL_STR = ("RPS", "RPL")
HB_STR = "^HB[^(P)]"

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)

def read_data(h5_file):
    logger.log(f"Reading data from {h5_file}", bcolors.WARNING)
    adata = sc.read_10x_h5(filename=h5_file)
    logger.log(f"Data from {h5_file} read", bcolors.OKBLUE)
    adata.var_names_make_unique()
    return adata

def prepare_4_scran(adata):
    adata_pp = adata.copy()
    sc.pp.normalize_total(adata_pp)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp, n_comps=15)
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(adata_pp, key_added="groups")

    data_mat = adata_pp.X.T
    # convert to CSC if possible. See https://github.com/MarioniLab/scran/issues/70
    if issparse(data_mat):
        if data_mat.nnz > 2**31 - 1:
            data_mat = data_mat.tocoo()
        else:
            data_mat = data_mat.tocsc()
    print(f"Data matrix converted to {type(data_mat)}")
    return data_mat, adata_pp.obs["groups"]

parser: argparse.ArgumentParser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, 
                                                          description='Merge filtered 10x data with Soupx preprocessing and with matrix without doublets') 

group1 = parser.add_argument_group('Input and output files')
group1.add_argument('-s', '--samplename', type = str, help='Sample name', default="")
group1.add_argument('-o', '--outdir', type = str, help='Output directory', default="")
group1.add_argument('-f', '--filtered', type = str, help='Filtered .h5 data from 10X', default="")
group1.add_argument('-x', '--soupx', type = str, help='SoupX .npz file', default="")
group1.add_argument('-d', '--doublet_score', type = str, help='.tsv file with doublet scores', default="")



def main():
    args = parser.parse_args()
    outdir: str = args.outdir
    filtered_h5: str = args.filtered
    samplename = args.samplename
    soupx_file = args.soupx
    doublet_score_file = args.doublet_score

    sample_dir = os.path.join(outdir, samplename)
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(sample_dir, exist_ok=True)

    logger.log(f"Processing sample: {samplename}", bcolors.OKBLUE)    
    adata = sc.read_h5ad(filtered_h5)

    soupX_counts = scipy.sparse.load_npz(soupx_file)
    doublet_score = pd.read_csv(doublet_score_file, sep="\t", header=0)
    print(doublet_score.head())

    adata.layers["counts"] = adata.X
    adata.layers["soupX_counts"] = soupX_counts.T
    adata.X = adata.layers["soupX_counts"]

    adata.obs["scDblFinder_score"] = doublet_score['doublet_score'].tolist()
    adata.obs["scDblFinder_class"] = doublet_score['doublet_class'].tolist()
 
    
    logger.log(f"All data merged", bcolors.OKBLUE)
    print(f"Total number of genes: {adata.n_vars}")
    print(f"Number of genes after cell filter: {adata.n_vars}")
    print(adata.obs.scDblFinder_class.value_counts())

    adata.write(os.path.join(sample_dir, f"{samplename}_soupX_MarkedDoublets.h5ad"))

    data_mat, groups = prepare_4_scran(adata)
    scipy.sparse.save_npz(os.path.join(sample_dir, f"{samplename}_prep_scran.npz"), data_mat)
    groups = groups.to_frame()
    groups.to_csv(os.path.join(sample_dir, f"{samplename}_scran_groups.csv"))


if __name__ == "__main__":
    main()