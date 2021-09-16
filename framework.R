# This file includes functions of the deconvolution benchmarking framework. 

source('benchmark.R')
source('deconvolution.R')

#' @title get_name
#' @description
#' This function generates a result file name by joining all the input parameters using dot ".". 
#' @details
#' all the parameters are joined together. Only the ".rds" suffix of the rds data file is removed. 
#' @param 
#' params is a list of characters. 
#' @return
#' This function returns a string, which can be used as a result file name. 
#' @example
#' name = get_name(args)
get_name<-function(params){
    n = ''
    for( i in params){
        j = unlist(strsplit(i, "/"))
        j = gsub('.rds','',dplyr::last(j))
        n = paste0(n, j, sep = ".")
    }
    return(n)
}

prepare_self<-function(param){

	dataset = param[1]
	number_cells = round(as.numeric(param[2]))
	
	if(param[3]=='T'){
        sampleCT = TRUE
    }else{
        sampleCT = FALSE
    }
    if(param[4]=='T'){
        propsample = TRUE
    }else{
        propsample = FALSE
    }
    if(param[5]=='2'){
        pbmode = 2
    }else{
        pbmode = 1
    }
   
	### Read single cell data and metadata
	X = read_data(dataset)
#     print('===0')
	training <- as.numeric(unlist(sapply(unique(colnames(X$data)), function(x) {
				sample(which(colnames(X$data) %in% x), X$cell_counts[x]/2) })))
	testing <- which(!1:ncol(X$data) %in% training)
#     print('===1')
	Xtrain = prepare_train(X$data[,training], X$original_cell_names[training])
#     print('===2')
	pDataC = X$pData[training,]
	test <- X$data[,testing]
	colnames(test) <- X$original_cell_names[testing]
	Xtest <- Generator(sce = test, phenoData = X$pData[testing,], sampleCT = sampleCT, propsample = propsample, Num.mixtures = 1000, pool.size = number_cells, mode = pbmode)
	return(list('Xtest'=Xtest, 'Xtrain'=Xtrain, 'pDataC'=pDataC))
}


prepare_cross<-function(param){

	dataset1 = param[1]
	dataset2 = param[2]
	number_cells = round(as.numeric(param[3]))
	
	if(param[4]=='T'){
        sampleCT = TRUE
    }else{
        sampleCT = FALSE
    }
    if(param[5]=='T'){
        propsample = TRUE
    }else{
        propsample = FALSE
    }
    if(param[6]=='2'){
        pbmode = 2
    }else{
        pbmode = 1
    }
	#-------------------------------------------------------
	### Read single cell data and metadata
	X1 = read_data(dataset1)
	X2 = read_data(dataset2)

	#-------------------------------------------------------
	### QC
	if(FALSE){
		X1<-QC(X1)
		X2<-QC(X2)
	}
	
	#-------------------------------------------------------
	### get intersection for the features
	to_keep = intersect(rownames(X1$data), rownames(X2$data))
	print(paste0('number of intersect features: ', length(to_keep)) )
	X1$data = X1$data[to_keep,]
	X2$data = X2$data[to_keep,]
	
	#-------------------------------------------------------
	### Prepare the reference data (on train data)
	Xtrain = prepare_train(X1$data, X1$original_cell_names)
	pDataC = X1$pData

	#-------------------------------------------------------
	### Generation of 1000 pseudo-bulk mixtures (T) (on test data)
	test <- X2$data
	colnames(test) <- X2$original_cell_names

	Xtest <- Generator(sce = test, phenoData = X2$pData, Num.mixtures = 10000, sampleCT = sampleCT, propsample = propsample, pool.size = number_cells, mode = pbmode)
	
	return(list('Xtest'=Xtest, 'Xtrain'=Xtrain, 'pDataC'=pDataC))
}

#' @title Framework
#' 
#' @description
#' This function contains the whole framework after data preprocessing and before deconvolution.
#' @details
#' This function is wrapped, because self_reference, cross_reference and bulk_reference can all reuse it. 
#' @param 
#' @param 
#' @param 
#' @param 
#' @param 
#' 
#' @return
#' The final deconvolution results. 
#' @example
#' 
Framework<-function(deconv_type,
					NormTrans,
					Xtest,
					Xtrain,
					normalization,
					normalization_scT,
					normalization_scC,
					transformation,
					marker_strategy,
					to_remove,
					P,
					method,
					pDataC){
	#-------------------------------------------------------
	### Transformation, scaling/normalization, marker selection for bulk deconvolution methods and deconvolution:
	if(deconv_type == "bulk"){
		if(NormTrans){
			T = Scaling(Xtest$T, normalization)
			C = Scaling(Xtrain$C, normalization)
			
			T = Transformation(T, transformation)
			C = Transformation(C, transformation)
		}else{
			T = Transformation(Xtest$T, transformation)
			C = Transformation(Xtrain$C, transformation)

			T = Scaling(T, normalization)
			C = Scaling(C, normalization)
		}

		# marker selection (on training data) 
		marker_distrib = marker_strategies(Xtrain$markers, marker_strategy, C)

		#If a cell type is removed, only meaningful mixtures where that CT was present (proportion < 0) are kept:
		if(to_remove != "none"){

			T <- T[,P[to_remove,] != 0]
			C <- C[, colnames(C) %in% rownames(P) & (!colnames(C) %in% to_remove)]
			P <- P[!rownames(P) %in% to_remove, colnames(T)]
			Xtrain$ref = Xtrain$ref[,colnames(Xtrain$ref) %in% rownames(P) & (!colnames(Xtrain$ref) %in% to_remove)]
			marker_distrib <- marker_distrib[marker_distrib$CT %in% rownames(P) & (marker_distrib$CT != to_remove),]
		}

	} else if (deconv_type == "sc"){
		if(NormTrans){
			T = Scaling(Xtest$T, normalization_scT)
			C = Scaling(Xtrain$train_cellID, normalization_scC)
			
			T = Transformation(T, transformation)
			C = Transformation(C, transformation)
		}else{
			T = Transformation(Xtest$T, transformation)
			C = Transformation(Xtrain$train_cellID, transformation)

			T = Scaling(T, normalization_scT)
			C = Scaling(C, normalization_scC)
		}
		#If a cell type is removed, only meaningful mixtures where that CT was present (proportion < 0) are kept:
		if(to_remove != "none"){

			T <- T[,P[to_remove,] != 0]
			C <- C[,pDataC$cellType != to_remove]
			P <- P[!rownames(P) %in% to_remove, colnames(T)]
			pDataC <- pDataC[pDataC$cellType != to_remove,]

		}
		marker_distrib = NULL
	}
	## 
	# T is the pseudo-bulk exprs
	# C is the reference exprs (average for each cell type)
	# P is the preassigned proportion for the pseudo-bulk
	# pDataC is the meta-data for the reference
#     saveRDS(T, 'T.rds')
#     saveRDS(C, 'C.rds')
#     saveRDS(pDataC, 'pDataC.rds')
	RESULTS = Deconvolution(T = T, 
							C = C, 
							method = method, 
							P = P, 
							elem = to_remove, 
							refProfiles.var = Xtrain$ref, 
							STRING = as.character(sample(1:10000, 1)), 
# 							STRING = "text", 
							marker_distrib = marker_distrib, 
							phenoDataC = pDataC) 
    
	return(RESULTS)
}

#' @title self_reference
#' 
#' @description
#' Use one single-cell references, split the dataset into two parts: train and test. 
#' The train is used as reference, while the test is used to make pseudobulk. 
#' @details
#' 
#' @param param: a list of 12 parameters:
#' 1: dataset name
#' 2: Transformation: none (defalt), log, sqrt, vst
#' 3: deconvolution type: bulk, sc
#' 4: Normalization for C, normalization
#' 5: Normalization for T, marker strategy
#' 6: Deconvolution method
#' 7: number of cells used
#' 8: remove cell type or not (none: default)
#' 9: number of cores used.
#' 10: sampleCT.
#' 11: propsample.
#' 12: Normalize first (T) or Transform first(F). 
#' 
#' @return
#' The deconvolution results as a table.
#' @example
#' 
self_reference<-function(param){
	
	if(length(param)!=13){

		print("Please check that all required parameters are indicated or are correct")
		print("Example usage for bulk deconvolution methods: 'Rscript Master_deconvolution.R baron none bulk TMM all nnls 100 none 1'")
		print("Example usage for single-cell deconvolution methods: 'Rscript Master_deconvolution.R baron none sc TMM TMM MuSiC 100 none 1'")
		stop()
	} 

	flag = FALSE

	dataset = param[1]
	transformation = param[2]
	deconv_type = param[3]

	if(deconv_type == "bulk"){
		normalization = param[4]
		marker_strategy = param[5]
	} else if (deconv_type == "sc") {
		normalization_scC = param[4]
		normalization_scT = param[5]
	} else {
		print("Please enter a valid deconvolution framework")
		stop()
	}

	method = param[6]
	number_cells = round(as.numeric(param[7]), digits = -2) #has to be multiple of 100
	to_remove = param[8]
	num_cores = min(as.numeric(param[9]),parallel::detectCores()-1)
	if(param[10]=='T'){
        sampleCT = TRUE
    }else{
        sampleCT = FALSE
    }
    if(param[11]=='T'){
        propsample = TRUE
    }else{
        propsample = FALSE
    }
    if(param[12]=='T'){
		NormTrans = TRUE
	}else{
        NormTrans = FALSE
    }
    
    ########### prepare_data ########################
	#-------------------------------------------------------
	### Read single cell data and metadata
	X = read_data(dataset)

	#-------------------------------------------------------
	### QC
	if(FALSE){
		X<-QC(X)
	}
	#-------------------------------------------------------
	### Data split into training/test  
	training <- as.numeric(unlist(sapply(unique(colnames(X$data)), function(x) {
				sample(which(colnames(X$data) %in% x), X$cell_counts[x]/2) })))
	testing <- which(!1:ncol(X$data) %in% training)
	
	#-------------------------------------------------------
	### Prepare the reference data (on train data)
	Xtrain = prepare_train(X$data[,training], X$original_cell_names[training])
	pDataC = X$pData[training,]

	#-------------------------------------------------------
	### Generation of 1000 pseudo-bulk mixtures (T) (on test data)
	test <- X$data[,testing]
	colnames(test) <- X$original_cell_names[testing]

	Xtest <- Generator(sce = test, phenoData = X$pData[testing,], sampleCT = sampleCT, propsample = propsample, Num.mixtures = 1000, pool.size = number_cells)
# 	P <- Xtest$P
	##################################################
	
	return(Framework(deconv_type,
					NormTrans,
					Xtest,
					Xtrain,
					normalization,
					normalization_scT,
					normalization_scC,
					transformation,
					marker_strategy,
					to_remove,
					Xtest$P,
					method,
					pDataC))
	
}

self_reference_pro<-function(param){
	
	if(length(param)!=12){

		print("Please check that all required parameters are indicated or are correct")
		print("Example usage for bulk deconvolution methods: 'Rscript Master_deconvolution.R baron none bulk TMM all nnls 100 none 1'")
		print("Example usage for single-cell deconvolution methods: 'Rscript Master_deconvolution.R baron none sc TMM TMM MuSiC 100 none 1'")
		stop()
	} 

	flag = FALSE

	dataset = param[1]
	transformation = param[2]
	deconv_type = param[3]

	if(deconv_type == "bulk"){
		normalization = param[4]
		marker_strategy = param[5]
	} else if (deconv_type == "sc") {
		normalization_scC = param[4]
		normalization_scT = param[5]
	} else {
		print("Please enter a valid deconvolution framework")
		stop()
	}

	method = param[6]
# 	number_cells = round(as.numeric(param[7]), digits = -2) #has to be multiple of 100
	to_remove = param[8]
	num_cores = min(as.numeric(param[9]),parallel::detectCores()-1)
# 	if(param[10]=='T'){
#         sampleCT = TRUE
#     }else{
#         sampleCT = FALSE
#     }
#     if(param[11]=='T'){
#         propsample = TRUE
#     }else{
#         propsample = FALSE
#     }
    if(param[12]=='T'){
		NormTrans = TRUE
	}else{
        NormTrans = FALSE
    }
    
    pbname = get_name(c(param[1], param[7], param[10], param[11], param[13]))
    
    x = read_bulk(sprintf('RDS/p.%s.rds', pbname))
	Xtest = x$Xtest
	pDataC = x$pDataC
	Xtrain = x$Xtrain
	
	return(Framework(deconv_type,
					NormTrans,
					Xtest,
					Xtrain,
					normalization,
					normalization_scT,
					normalization_scC,
					transformation,
					marker_strategy,
					to_remove,
					Xtest$P,
					method,
					pDataC))
	
}

#' @title cross_reference
#' 
#' @description
#' Use two single-cell references, one to build pseudobulk, the other one to use as deconvolution reference. 
#' @details
#' 
#' @param param: a list of 13 parameters:
#' 1: dataset1 name for reference
#' 2: dataset2 name for pseudobulk
#' 3: Transformation: none (defalt), log, sqrt, vst
#' 4: deconvolution type: bulk, sc
#' 5: Normalization for C, normalization
#' 6: Normalization for T, marker strategy
#' 7: Deconvolution method
#' 8: number of cells used
#' 9: remove cell type or not (none: default)
#' 10: number of cores used.
#' 11: sampleCT.
#' 12: propsample.
#' 13: Normalize first (T) or Transform first(F). 
#' 
#' @return
#' The deconvolution results as a table.
#' @example
#' 
cross_reference<-function(param){
	
	if(length(param)!=13){

		print("Please check that all required parameters are indicated or are correct")
		print("Example usage for bulk deconvolution methods: 'Rscript Master_deconvolution.R baron none bulk TMM all nnls 100 none 1'")
		print("Example usage for single-cell deconvolution methods: 'Rscript Master_deconvolution.R baron none sc TMM TMM MuSiC 100 none 1'")
		stop()
	} 

	flag = FALSE

	dataset1 = param[1]
	dataset2 = param[2]
	transformation = param[3]
	deconv_type = param[4]

	if(deconv_type == "bulk"){
		normalization = param[5]
		marker_strategy = param[6]
	} else if (deconv_type == "sc") {
		normalization_scC = param[5]
		normalization_scT = param[6]
	} else {
		print("Please enter a valid deconvolution framework")
		stop()
	}

	method = param[7]
	number_cells = round(as.numeric(param[8]), digits = -2) #has to be multiple of 100
	to_remove = param[9]
	num_cores = min(as.numeric(param[10]),parallel::detectCores()-1)
    if(param[11]=='T'){
        sampleCT = TRUE
    }else{
        sampleCT = FALSE
    }
    if(param[12]=='T'){
        propsample = TRUE
    }else{
        propsample = FALSE
    }
    if(param[13]=='T'){
		NormTrans = TRUE
	}else{
        NormTrans = FALSE
    }

	#-------------------------------------------------------
	### Read single cell data and metadata
	X1 = read_data(dataset1)
	X2 = read_data(dataset2)

	#-------------------------------------------------------
	### QC
	if(FALSE){
		X1<-QC(X1)
		X2<-QC(X2)
	}
	
	#-------------------------------------------------------
	### get intersection for the features
	to_keep = intersect(rownames(X1$data), rownames(X2$data))
	print(paste0('number of intersect features: ', length(to_keep)) )
	X1$data = X1$data[to_keep,]
	X2$data = X2$data[to_keep,]
	
	#-------------------------------------------------------
	### Prepare the reference data (on train data)
	Xtrain = prepare_train(X1$data, X1$original_cell_names)
	pDataC = X1$pData

	#-------------------------------------------------------
	### Generation of 1000 pseudo-bulk mixtures (T) (on test data)
	test <- X2$data
	colnames(test) <- X2$original_cell_names

	Xtest <- Generator(sce = test, phenoData = X2$pData, Num.mixtures = 1000, sampleCT = sampleCT, propsample = propsample, pool.size = number_cells)
# 	P <- Xtest$P

	return(Framework(deconv_type,
					NormTrans,
					Xtest,
					Xtrain,
					normalization,
					normalization_scT,
					normalization_scC,
					transformation,
					marker_strategy,
					to_remove,
					Xtest$P,
					method,
					pDataC))
}

cross_reference_pro<-function(param){
	
	if(length(param)!=13){

		print("Please check that all required parameters are indicated or are correct")
		print("Example usage for bulk deconvolution methods: 'Rscript Master_deconvolution.R baron none bulk TMM all nnls 100 none 1'")
		print("Example usage for single-cell deconvolution methods: 'Rscript Master_deconvolution.R baron none sc TMM TMM MuSiC 100 none 1'")
		stop()
	} 

	flag = FALSE

	dataset1 = param[1]
	dataset2 = param[2]
	transformation = param[3]
	deconv_type = param[4]

	if(deconv_type == "bulk"){
		normalization = param[5]
		marker_strategy = param[6]
	} else if (deconv_type == "sc") {
		normalization_scC = param[5]
		normalization_scT = param[6]
	} else {
		print("Please enter a valid deconvolution framework")
		stop()
	}

	method = param[7]
# 	number_cells = round(as.numeric(param[8]), digits = -2) #has to be multiple of 100
	to_remove = param[9]
	num_cores = min(as.numeric(param[10]),parallel::detectCores()-1)
#     if(param[11]=='T'){
#         sampleCT = TRUE
#     }else{
#         sampleCT = FALSE
#     }
#     if(param[12]=='T'){
#         propsample = TRUE
#     }else{
#         propsample = FALSE
#     }
    if(param[13]=='T'){
		NormTrans = TRUE
	}else{
        NormTrans = FALSE
    }

	pbname = get_name(c(param[1], param[2], param[8], param[11], param[12], param[14]))
    
    x = read_bulk(sprintf('RDS/q.%s.rds', pbname))
	Xtest = x$Xtest
	pDataC = x$pDataC
	Xtrain = x$Xtrain

	return(Framework(deconv_type,
					NormTrans,
					Xtest,
					Xtrain,
					normalization,
					normalization_scT,
					normalization_scC,
					transformation,
					marker_strategy,
					to_remove,
					Xtest$P,
					method,
					pDataC))
}

#' @title bulk_reference
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
bulk_reference<-function(param){
	
	if(length(param)!=11){

		print("Please check that all required parameters are indicated or are correct")
		print("Example usage for bulk deconvolution methods: 'Rscript Master_deconvolution.R baron none bulk TMM all nnls 100 none 1'")
		print("Example usage for single-cell deconvolution methods: 'Rscript Master_deconvolution.R baron none sc TMM TMM MuSiC 100 none 1'")
		stop()
	} 

	flag = FALSE

	bulk = param[1]
	dataset = param[2]
	transformation = param[3]
	deconv_type = param[4]

	if(deconv_type == "bulk"){
		normalization = param[5]
		marker_strategy = param[6]
	} else if (deconv_type == "sc") {
		normalization_scC = param[5]
		normalization_scT = param[6]
	} else {
		print("Please enter a valid deconvolution framework")
		stop()
	}

	method = param[7]
	number_cells = round(as.numeric(param[8]), digits = -2) #has to be multiple of 100
	to_remove = param[9]
	num_cores = min(as.numeric(param[10]),parallel::detectCores()-1)
    if(param[11]=='T'){
		NormTrans = TRUE
	}else{
        NormTrans = FALSE
    }

	#-------------------------------------------------------
	### Read single cell data and metadata
	X1 = read_bulk(bulk)
	X2 = read_data(dataset)

	#-------------------------------------------------------
	### QC
	if(FALSE){
#  		X1<-QC(X1)
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

	#-------------------------------------------------------
	### Generation of 1000 pseudo-bulk mixtures (T) (on test data)
#  	test <- X2$data
#  	colnames(test) <- X2$original_cell_names

#  	Xtest <- Generator(sce = test, phenoData = X2$pData, Num.mixtures = 1000, pool.size = number_cells)
#  	P <- Xtest$P

	P<-as.data.frame(matrix(0.1, nrow=length(unique(colnames(X2$data))), ncol=dim(X1$data)[2]))
	rownames(P)<-unique(colnames(X2$data))
	colnames(P)<-colnames(X1$data)
	Xtest<-list('T'=X1$data, 'P'=P)

	return(Framework(deconv_type,
					NormTrans,
					Xtest,
					Xtrain,
					normalization,
					normalization_scT,
					normalization_scC,
					transformation,
					marker_strategy,
					to_remove,
					Xtest$P,
					method,
					pDataC))
}

#' @title bulk_2references
#' 
#' @description
#' Use TWO single-cell references to deconvolve a bulk data. And compare the deconvolution results. 
#' @details
#' 
#' @param param: a list of 12 parameters:
#' 1: bulk
#' 2: reference1
#' 3: reference2
#' 4: Transformation: none (defalt), log, sqrt, vst
#' 5: deconvolution type: bulk, sc
#' 6: Normalization for C, normalization
#' 7: Normalization for T, marker strategy
#' 8: Deconvolution method
#' 9: number of cells used
#' 10: remove cell type or not (none: default)
#' 11: number of cores used.
#' 12: Normalize first (T) or Transform first(F).
#' 
#' @return
#' The combined deconvolution results as a table.
#' @example
#' 
bulk_2references<-function(param){
	
	if(length(param)!=11){

		print("Please check that all required parameters are indicated or are correct")
		print("Example usage for bulk deconvolution methods: 'Rscript Master_deconvolution.R baron none bulk TMM all nnls 100 none 1'")
		print("Example usage for single-cell deconvolution methods: 'Rscript Master_deconvolution.R baron none sc TMM TMM MuSiC 100 none 1'")
		stop()
	} 

	flag = FALSE

	### arguments

	
	param1<-c(param[1], param[2], param[4], param[5], param[6], param[7], param[8], param[9], param[10], param[11], param[12])
	param2<-c(param[1], param[3], param[4], param[5], param[6], param[7], param[8], param[9], param[10], param[11], param[12])
	
	RESULTS1 = bulk_reference(param1)
	RESULTS2 = bulk_reference(param2)
    
    RESULTS = merge(RESULTS1, RESULTS2, by=c("CT","tissue"))

	return(RESULTS)
}

combine_reference_pro<-function(param){
	
	if(length(param)!=9){
#         print(length(param))
#         print(param)
		print("Please check that all required parameters are indicated or are correct")
		print("Example usage for bulk deconvolution methods: 'Rscript Master_deconvolution.R baron none bulk TMM all nnls 100 none 1'")
		print("Example usage for single-cell deconvolution methods: 'Rscript Master_deconvolution.R baron none sc TMM TMM MuSiC 100 none 1'")
		stop()
	} 

	dataset = param[1]
	transformation = param[2]
	deconv_type = param[3]

	if(deconv_type == "bulk"){
		normalization = param[4]
		marker_strategy = param[5]
	} else if (deconv_type == "sc") {
		normalization_scC = param[4]
		normalization_scT = param[5]
	} else {
		print("Please enter a valid deconvolution framework")
		stop()
	}

	method = param[6]
	to_remove = param[7]
	num_cores = min(as.numeric(param[8]),parallel::detectCores()-1)

    if(param[9]=='T'){
		NormTrans = TRUE
	}else{
        NormTrans = FALSE
    }
    
#     x = read_bulk(sprintf('RDS/%s.rds', param[1]))
    x = read_bulk(param[1])
	Xtest = x$Xtest
	pDataC = x$pDataC
	Xtrain = x$Xtrain

	return(Framework(deconv_type,
					NormTrans,
					Xtest,
					Xtrain,
					normalization,
					normalization_scT,
					normalization_scC,
					transformation,
					marker_strategy,
					to_remove,
					Xtest$P,
					method,
					pDataC))
}
