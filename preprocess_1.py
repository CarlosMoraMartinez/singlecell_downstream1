# %%
import argparse
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import seaborn as sns
from scipy.stats import median_abs_deviation
import os
from myutils import *
import scipy.sparse
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

def calculate_qc_metrics(adata):
    adata.var["mt"] = adata.var_names.str.startswith(MIT_STR)
    adata.var["ribo"] = adata.var_names.str.startswith(RIBOSOMAL_STR)
    adata.var["hb"] = adata.var_names.str.contains(HB_STR)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )
    return adata

def plot_qc_metrics(adata, outdir, name):
    p1 = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
    p1.savefig(os.path.join(outdir, f"{name}_total_counts.png"))
    plt.close(p1.fig)

    p2 = sc.pl.violin(adata, "pct_counts_mt", show=False)
    p2.figure.savefig(os.path.join(outdir, f"{name}_pct_counts_mt.png"))
    plt.close(p2.figure)

    p3 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False)
    p3.figure.savefig(os.path.join(outdir, f"{name}_total_counts_vs_n_genes_by_counts.png"))
    plt.close(p3.figure)

def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

def filter_cells(adata, outlier_median_desv_coverages = 5, outlier_median_desv_mitochnondrial = 8):
    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", outlier_median_desv_coverages)
        | is_outlier(adata, "log1p_n_genes_by_counts", outlier_median_desv_coverages)
        | is_outlier(adata, "pct_counts_in_top_20_genes", outlier_median_desv_coverages)
    )
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
        adata.obs["pct_counts_mt"] > outlier_median_desv_mitochnondrial
    )
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
    return adata

def preprocess_for_soupx(adata):
    adata_pp = adata.copy()
    sc.pp.normalize_per_cell(adata_pp)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp)
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(adata_pp, key_added="soupx_groups")
    soupx_groups = adata_pp.obs["soupx_groups"]
    return adata_pp, soupx_groups

parser: argparse.ArgumentParser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, 
                                                          description='Preprocess 10x data for Soupx')

group1 = parser.add_argument_group('Input and output files')
group1.add_argument('-s', '--samplename', type = str, help='Sample name', default="")
group1.add_argument('-o', '--outdir', type = str, help='Output directory', default="")
group1.add_argument('-f', '--filtered', type = str, help='Filtered .h5 data from 10X', default="")
group1.add_argument('-r', '--raw', type = str, help='Raw .h5 data from 10X', default="")

group2 = parser.add_argument_group('Input and output files')
group2.add_argument('-c', '--median_desv_coverages', type = int, help='Number of median abs deviations to detect outliers based on total_counts, n_genes_by_counts and pct_counts_in_top_20_genes', default=5)
group2.add_argument('-m', '--median_desv_mitochnondrial', type = int, help='Number of median abs deviations to detect outliers based on mitochondrial trancripts content', default=8)


def main():
    args = parser.parse_args()
    outdir: str = args.outdir
    raw_h5: str = args.raw
    filtered_h5: str = args.filtered
    samplename = args.samplename

    outlier_median_desv_coverages: int = args.median_desv_coverages
    outlier_median_desv_mitochnondrial: int = args.median_desv_mitochnondrial
    sample_dir = os.path.join(outdir, samplename)
    
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(sample_dir, exist_ok=True)

    logger.log(f"Processing sample: {samplename}", bcolors.OKBLUE)    
    adata = read_data(filtered_h5)
    adata_raw = read_data(raw_h5)

    adata = calculate_qc_metrics(adata)
    plot_qc_metrics(adata, sample_dir,  "raw")
    logger.log(f"Total number of cells: {adata.n_obs}", bcolors.OKBLUE)

    adata = filter_cells(adata, outlier_median_desv_coverages, outlier_median_desv_mitochnondrial)
    plot_qc_metrics(adata, sample_dir, "filtered")
    logger.log(f"Number of cells after filtering of low quality cells: {adata.n_obs}", bcolors.OKBLUE)

    adata_pp, soupx_groups = preprocess_for_soupx(adata)
    ## Save the filtered data
    adata.write_h5ad(os.path.join(sample_dir, f"{samplename}_filt.h5ad"))

    ## Save preprocessed data
    adata_pp.write_h5ad(os.path.join(sample_dir, f"{samplename}_preprocess.h5ad"))
    logger.log(f"Filtered data saved in {filtered_h5}", bcolors.OKGREEN)

    # write the output files for soupx

    data_tod = adata_raw.X.T
    scipy.sparse.save_npz(os.path.join(sample_dir, f"{samplename}_rawdata.npz"), data_tod)

    cells = pd.DataFrame(adata.obs_names, columns=["cell"])
    cells.to_csv(os.path.join(sample_dir, f"{samplename}_cells.csv"), index=False)

    genes = pd.DataFrame(adata.var_names, columns=["gene"])
    genes.to_csv(os.path.join(sample_dir, f"{samplename}_genes.csv"), index=False)

    data = adata.X.T
    scipy.sparse.save_npz(os.path.join(sample_dir, f"{samplename}_filtdata.npz"), data)

    soupx_groups = soupx_groups.to_frame()
    soupx_groups.to_csv(os.path.join(sample_dir, f"{samplename}_soupx_groups.csv"))

    
    logger.log(f"All data exported", bcolors.OKBLUE)



if __name__ == "__main__":
    main()

