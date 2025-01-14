import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import os
from myutils import *
import numpy as np
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

def normalize_log1p(adata):
    scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
    # log1p transform
    adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)
    return adata

def normalize_pearsonres(adata):
    analytic_pearson = sc.experimental.pp.normalize_pearson_residuals(adata, inplace=False)
    adata.layers["analytic_pearson_residuals"] = csr_matrix(analytic_pearson["X"])
    return adata

def normalize_scran(adata, scran_sf):
    adata.obs["size_factors"] = scran_sf["size_factors"].tolist()
    scran = adata.X / adata.obs["size_factors"].values[:, None]  
    adata.layers["scran_normalization"] = csr_matrix(sc.pp.log1p(scran))
    return adata

def plot_norm_counts(adata, outdir, name, object_name="log1p_norm", title_name="Shifted logarithm", bins=100):
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    p1 = sns.histplot(adata.obs["total_counts"], bins=bins, kde=False, ax=axes[0])
    axes[0].set_title("Total counts")
    p2 = sns.histplot(adata.layers[object_name].sum(1), bins=bins, kde=False, ax=axes[1])
    axes[1].set_title(title_name)
    fig.savefig(os.path.join(outdir, f"{name}_count_distribution_{object_name}.png"))
    plt.close(fig)

def get_deviance(adata, binomial_deviance, num_genes=4000):
    idx = binomial_deviance.argsort()[-num_genes:]
    mask = np.zeros(adata.var_names.shape, dtype=bool)
    mask[idx] = True
    
    adata.var["highly_deviant"] = mask
    adata.var["binomial_deviance"] = binomial_deviance
    return adata

def plot_highly_deviant_genes(adata, outdir, name):
    sc.pp.highly_variable_genes(adata, layer="scran_normalization")
    ax = sns.scatterplot(
    data=adata.var, x="means", y="dispersions", hue="highly_deviant", s=5
    )
    ax.set_xlim(None, 1.5)
    ax.set_ylim(None, 3)
    ax.figure.savefig(os.path.join(outdir, f"{name}_highly_deviant_genes.png"))
    plt.close('all')

def make_pca(adata, outdir, name):
    adata.var["highly_variable"] = adata.var["highly_deviant"]
    sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)
    plt.figure()
    #fig, axes = plt.subplots()
    p = sc.pl.pca_scatter(adata, color="total_counts", save=f"{name}_highly_deviant_genes")
    #p.savefig(os.path.join(outdir, f"{name}_highly_deviant_genes_PCA.png"))
    #print(type(p))
    plt.close('all')

parser: argparse.ArgumentParser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, 
                                                          description='Merge filtered 10x data with Soupx preprocessing and with matrix without doublets') 

group1 = parser.add_argument_group('Input and output files')
group1.add_argument('-s', '--samplename', type = str, help='Sample name', default="")
group1.add_argument('-o', '--outdir', type = str, help='Output directory', default="")
group1.add_argument('-f', '--filtered', type = str, help='Filtered .h5 data from 10X', default="")
group1.add_argument('-d', '--scran_sf', type = str, help='.tsv file with scran size factors', default="")
group1.add_argument('-c', '--deviance', type = str, help='.tsv file with scry deviance', default="")
group1.add_argument('-n', '--num_genes', type = int, help='Number of genes to consider as highly deviant', default=4000)

def main():
    args = parser.parse_args()
    outdir: str = args.outdir
    filtered_h5: str = args.filtered
    samplename = args.samplename
    scran_sf_file = args.scran_sf
    deviance_file = args.deviance
    num_genes = args.num_genes

    sample_dir = os.path.join(outdir, samplename)
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(sample_dir, exist_ok=True)

    logger.log(f"Processing sample: {samplename}", bcolors.OKBLUE)    
    adata = sc.read_h5ad(filtered_h5)

    scran_size_factors = pd.read_csv(scran_sf_file, sep="\t", header=0)
    #print(scran_size_factors.head())

    logger.log(f"Shifted log normalization", bcolors.OKBLUE)    
    adata = normalize_log1p(adata)
    plot_norm_counts(adata, sample_dir, samplename)
    
    # No memory
    #logger.log(f"Analytic Pearson residuals normalization", bcolors.OKBLUE)    
    #adata = normalize_pearsonres(adata)
    #plot_norm_counts(adata, sample_dir, samplename, "analytic_pearson_residuals", "Analytic Pearson residuals")

    logger.log(f"Adding Scran normalization", bcolors.OKBLUE)    
    adata = normalize_scran(adata, scran_size_factors)
    plot_norm_counts(adata, sample_dir, samplename, "scran_normalization", "log1p with Scran estimated size factors")

    deviance = pd.read_csv(deviance_file, sep="\t", header=0)
    adata = get_deviance(adata, deviance["binomial_deviance"].values, num_genes)

    adata.write(os.path.join(sample_dir, f"{samplename}_soupX_MarkedDoublets_Norm.h5ad"))

    plot_highly_deviant_genes(adata, sample_dir, f"{samplename}_{num_genes}")
    #make_pca(adata, sample_dir, f"{samplename}_{num_genes}")
    



if __name__ == "__main__":
    main()