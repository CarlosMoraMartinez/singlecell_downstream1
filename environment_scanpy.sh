

mamba create -n singlecell1 python=3.12
mamba activate singlecell1
mamba install -y -c conda-forge numpy 
mamba install -y -c conda-forge pandas 
conda install -y conda-forge::matplotlib
mamba install -y -c conda-forge seaborn
mamba install -y -c conda-forge scipy
mamba install -y -c conda-forge scikit-learn
mamba install -y -c conda-forge anndata 
mamba install -y -c conda-forge scanpy 
mamba install -y -c conda-forge mudata
mamba install -y -c conda-forge jupyterlab
conda install -c conda-forge python-igraph

mamba install -y -c conda-forge rpy2.robjects 
mamba install -y -c conda-forge rpy2
mamba install -y -c conda-forge anndata2ri