# This file includes functions of the deconvolution pipeline. 

source('benchmark.R')
source('deconvolution.R')
source('framework.R')

#' @title pipeline
#' 
#' @description
#' Use a single-cell references to deconvolve a bulk data. 
#' @details
#' 
#' @param param: a list of 11 parameters:
#' 1: bulk name
#' 2: reference dataset name
#' 3: Transformation: none (defalt), log, sqrt, vst
#' 4: deconvolution type: bulk, sc
#' 5: Normalization for C, normalization
#' 6: Normalization for T, marker strategy
#' 7: Deconvolution method
#' 8: number of cells used
#' 9: remove cell type or not (none: default)
#' 10: number of cores used.
#' 11: Normalize first (T) or Transform first(F).
#' 
#' @return
#' The deconvolution results as a table.
#' @example
#' 
pipeline<-function(param){
	
	if(length(param)!=11){

		print("Please check that all required parameters are indicated or are correct")
		print("Example usage for bulk deconvolution methods: 'Rscript Master_deconvolution.R baron none bulk TMM all nnls 100 none 1'")
		print("Example usage for single-cell deconvolution methods: 'Rscript Master_deconvolution.R baron none sc TMM TMM MuSiC 100 none 1'")
		stop()
	} 

	flag = FALSE

	bulk = param[1]
	dataset = param[2]

	method = param[3]
	
	bulk_methods = c("CIBERSORT","DeconRNASeq","OLS","nnls","FARDEEP","RLR","DCQ","elasticNet","lasso","ridge","EPIC",
					 "DSA","ssKL","ssFrobenius","dtangle", "deconf", "proportionsInAdmixture", "svmdecon", "EpiDISH","CAMmarker","CDSeq" )
	sc_methods = c("MuSiC","BisqueRNA","DWLS","deconvSeq","SCDC","bseqsc","CPM","TIMER")

	if(method %in% bulk_methods){
		deconv_type = "bulk"
		normalization = param[4]
		marker_strategy = 'all'
	} else if (method %in% sc_methods) {
		deconv_type = "sc"
		normalization_scC = param[4]
		normalization_scT = param[4]
	} else {
		print("Please enter a valid deconvolution framework")
		stop()
	}

	

	#-------------------------------------------------------
	### Read single cell data and metadata
	X1 = read_bulk(bulk)
	X2 = read_data(dataset)

	#-------------------------------------------------------
	### QC
	if(FALSE){
		X2<-QC(X2)
	}
	
	#-------------------------------------------------------
	### get intersection for the features
	to_keep = intersect(rownames(X1$data), rownames(X2$data))
	print(paste0('number of intersect features: ', length(to_keep)) )
	X1$data = X1$data[to_keep,]
	X2$data = X2$data[to_keep,]
	
	#-------------------------------------------------------
	### Prepare the reference data
	Xtrain = prepare_train(X2$data, X2$original_cell_names)
	pDataC = X2$pData

	P<-as.data.frame(matrix(0.1, nrow=length(unique(colnames(X2$data))), ncol=dim(X1$data)[2]))
	rownames(P)<-unique(colnames(X2$data))
	colnames(P)<-colnames(X1$data)
	Xtest<-list('T'=X1$data, 'P'=P)

	return(Framework(deconv_type,
					NormTrans = TRUE,
					Xtest,
					Xtrain,
					normalization,
					normalization_scT,
					normalization_scC,
					transformation = 'none',
					marker_strategy,
					to_remove = 'none',
					Xtest$P,
					method,
					pDataC))
}
