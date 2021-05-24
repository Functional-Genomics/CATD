args <- commandArgs(trailingOnly=TRUE)

# source('helper_functions.R')
# anndata2expressionset(args[1], args[2])


Sys.setenv(MKL_SERVICE_FORCE_INTEL=1)
library(sceasy)
library(reticulate)
use_condaenv('r4-base')
use_python("/nfs/leia/research/ma/chichau/bin/anaconda3/bin/python", required = TRUE)
# py_config()
loompy <- reticulate::import('loompy')

sceasy::convertFormat(args[1], from="anndata", to="seurat",
                       outFile=args[2])
