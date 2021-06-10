# This file includes functions to convert the file format from python-based AnnData (https://anndata.readthedocs.io/en/latest/anndata.AnnData.html)
# to R-based Seurat (https://satijalab.org/seurat/) format or ExpressionSet (https://www.rdocumentation.org/packages/Biobase/versions/2.32.0/topics/ExpressionSet) format. 

#' @title read_h5
#' @description
#' This function reads the AnnData h5ad file and save the .X, .obs and .raw.X data as a list. 
#' @details
#' var_names are saved as rownames, obs_names are saved as colnames.
#' @param 
#' h5file: the AnnData h5ad file. 
#' 
#' @return
#' A list of 3 values: X is the .X of AnnData, rawX is the .raw.X of AnnData, and obs is the .obs of AnnData. 
#' @example
#' x <- read_h5("example.h5")
read_h5<-function(h5file){
    suppressPackageStartupMessages(library(reticulate))
    suppressPackageStartupMessages(library(Matrix))
    handle <- import("anndata", convert = FALSE)
    ad <- handle$read_h5ad(h5file)

    obs <- py_to_r(ad$obs)
    if ("n_counts" %in% colnames(x = obs)) {
    colnames(x = obs) <- gsub(
          pattern = "n_counts",
          replacement = "nUMI",
          x = colnames(x = obs)
        )
    }

    genenames <- rownames(py_to_r(ad$var))
    cellnames <- rownames(obs)

    ncells = length(cellnames)
    ngenes = length(genenames)
    # Read X
    if("indptr" %in% names(ad$X)){
        ad.X <- sparseMatrix(
            i = as.numeric(x = ad$X$indices),
            p = as.numeric(x = ad$X$indptr),
            x = as.numeric(x = ad$X$data),
            index1 = FALSE,
            dims=c(ngenes,
                   ncells)
        )
    }else{
        ad.X <- t(py_to_r(ad$X))
    }
    rownames(ad.X) <- genenames
    colnames(ad.X) <- cellnames

    # Read raw
    ad.raw.X = NA
    if("raw" %in% names(ad)){
        if("X" %in% names(ad$raw)){
            genenames <- rownames(py_to_r(ad$raw$var))
            ngenes = length(genenames)
            if("indptr" %in% names(ad$raw$X)){
                ad.raw.X <- sparseMatrix(
                    i = as.numeric(x = ad$raw$X$indices),
                    p = as.numeric(x = ad$raw$X$indptr),
                    x = as.numeric(x = ad$raw$X$data),
                    index1 = FALSE,
                    dims=c(ngenes,
                           ncells)
                )
            }else{
                ad.raw.X <- t(py_to_r(ad$raw$X))
            }
            rownames(ad.raw.X) <- genenames
            colnames(ad.raw.X) <- cellnames
        }else{
            genenames <- rownames(py_to_r(ad$var))
            ngenes = length(genenames)
            if("indptr" %in% names(ad$X)){
                ad.raw.X <- sparseMatrix(
                    i = as.numeric(x = ad$X$indices),
                    p = as.numeric(x = ad$X$indptr),
                    x = as.numeric(x = ad$X$data),
                    index1 = FALSE,
                    dims=c(ngenes,
                           ncells)
                )
            }else{
                ad.raw.X <- t(py_to_r(ad$X))
            }
            rownames(ad.raw.X) <- genenames
            colnames(ad.raw.X) <- cellnames
        }
    }
    return(list(X=ad.X, rawX=ad.raw.X, obs=obs))
}

#' @title anndata2expressionset
#' 
#' @description
#' Convert the AnnData h5ad file to ExpressionSet file.
#' @details
#' 
#' @param h5file: the AnnData h5ad file.
#' @param esfile: the ExpressionSet file to be saved.
#' 
#' @return
#' NA
#' @example
#' anndata2expressionset('anndata.h5', 'expr.rds')
anndata2expressionset<-function(h5file, esfile){
	require(Biobase)
    df = read_h5(h5file)
    ad = Biobase::ExpressionSet(assayData = as.matrix(df$X),
								phenoData = Biobase::AnnotatedDataFrame(df$obs))
    saveRDS(ad, esfile)
}

# main
# usage: !/nfs/leia/research/ma/chichau/bin/anaconda3/envs/r4-base/bin/Rscript format.R input.h5 output.rds
# use sceasy(https://github.com/cellgeni/sceasy) to convert the Anndata file to Seurat file. 
args <- commandArgs(trailingOnly=TRUE)
library(sceasy)
library(reticulate)
loompy <- reticulate::import('loompy')
sceasy::convertFormat(args[1], 
					  from="anndata", 
					  to="seurat",
                      outFile=args[2])
