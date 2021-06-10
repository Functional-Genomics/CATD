# install latest anaconda

conda update conda
conda update anaconda

# TOC
conda install -c conda-forge nodejs
pip install jupyter_contrib_nbextensions
jupyter labextension install @jupyterlab/toc

pip install rpy2

# install scvi-tools
pip install scvi-tools

pip install sccaf

conda install -c conda-forge cmake
conda install -c conda-forge r-biocmanager

pip install MulticoreTSNE

conda install -c conda-forge r-base



conda create --name r4-base
conda activate r4-base
conda install -c conda-forge r-base
conda install -c conda-forge/label/gcc7 r-base

conda install -c conda-forge r-seuratdisk

conda install -c conda-forge r-irkernel


R
install.packages('IRkernel')
IRkernel::installspec()  # to register the kernel in the current R installation
jupyter labextension install @techrah/text-shortcuts  # for RStudio’s shortcuts


conda create --name firstEnv
conda install -c anaconda ipykernel
python -m ipykernel install --user --name=firstEnv


conda install -c bioconda r-sceasy

gene programming
