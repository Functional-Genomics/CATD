
#############################
# Conversion

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

anndata2expressionset<-function(h5file, save){
    df = read_h5(h5file)
    ad = Biobase::ExpressionSet(assayData = as.matrix(df$X),phenoData = Biobase::AnnotatedDataFrame(df$obs))
    saveRDS(ad, save)
}
