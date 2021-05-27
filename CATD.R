#-------------------------------------------------------
### Helper functions + CIBERSORT external code
source('./helper_functions.R')
set.seed(24)
require(limma); require(dplyr); require(pheatmap); require(Matrix)
args <- commandArgs(trailingOnly=TRUE)

read_data<-function(dataset){
    
    ad = readRDS(dataset)
    if (class(ad)== "Seurat") {
#        print("input of class SeuratObject")
        data = ad@assays$RNA@counts 
        pData = ad@meta.data
    } else if (class(ad) =="ExpressionSet"){
        
#        print("input of class ExpressionSet")
        data = ad@assayData$exprs
        pData = ad@phenoData@data
               
    } else {
		print("the format of the input data is neither Seurat nor ExpressionSet")
    }
    rownames(pData)<-NULL

	original_cell_names = colnames(data)
	colnames(data) <- as.character(pData$cellType[match(colnames(data),pData$cellID)])
	cell_counts = table(colnames(data))
	return(list('data'=data, 'pData'=pData, 'original_cell_names'=original_cell_names, 'cell_counts'=cell_counts))
}

read_bulk<-function(dataset){
	data = readRDS(dataset)
	return(list('data'=data))
}

QC_feature<-function(data){
	
	keep <- which(Matrix::rowSums(data > 0) >= 3)
	data = data[keep,]
}

QC<-function(X){
	# filter features
	keep <- which(Matrix::rowSums(X$data > 0) >= 3)
	X$data = X$data[keep,]
	# Keep CTs with >= 4 cells after QC
	to_keep = names(X$cell_counts)[X$cell_counts >= 4]
	X$pData <- X$pData[X$pData$cellType %in% to_keep,]
	to_keep = which(colnames(X$data) %in% to_keep)   
	X$data <- X$data[,to_keep]
	X$original_cell_names <- X$original_cell_names[to_keep]
	return(X)
}

ref_mat<-function(train){
	cellType <- colnames(train)
	group = list()
	for(i in unique(cellType)){ 
		group[[i]] <- which(cellType %in% i)
	}
	#C should be made with the mean (not sum) to agree with the way markers were selected
	C = lapply(group,function(x) Matrix::rowMeans(train[,x])) 
	C = round(do.call(cbind.data.frame, C))

	refProfiles.var = lapply(group,function(x) train[,x])
	refProfiles.var = lapply(refProfiles.var, function(x) matrixStats::rowSds(Matrix::as.matrix(x)))
	refProfiles.var = round(do.call(cbind.data.frame, refProfiles.var))
	rownames(refProfiles.var) <- rownames(train)
	return(list('C'=C, 'ref'=refProfiles.var))
}


marker_selection<-function(train){
	#-------------------------------------------------------
	#Normalization of "train" followed by marker selection 
	#for marker selection, keep genes where at least 30% of cells within a cell type have a read/UMI count different from 0
	cellType = colnames(train) 
	keep <- sapply(unique(cellType), function(x) {
		CT_hits = which(cellType %in% x)
		size = ceiling(0.3*length(CT_hits))
		Matrix::rowSums(train[,CT_hits,drop=FALSE] != 0) >= size
	})
	train = train[Matrix::rowSums(keep) > 0,]
	train2 = Normalization(train)

	# INITIAL CONTRASTS for marker selection WITHOUT taking correlated CT into account 
	#[compare one group with average expression of all other groups]
	annotation = factor(colnames(train2))
	design <- model.matrix(~0+annotation)
	colnames(design) <- unlist(lapply(strsplit(colnames(design),"annotation"), function(x) x[2]))
	cont.matrix <- matrix((-1/ncol(design)),nrow=ncol(design),ncol=ncol(design))
	colnames(cont.matrix) <- colnames(design)
	diag(cont.matrix) <- (ncol(design)-1)/ncol(design)

	v <- limma::voom(train2, design=design, plot=FALSE) 
	fit <- limma::lmFit(v, design)
	fit2 <- limma::contrasts.fit(fit, cont.matrix)
	fit2 <- limma::eBayes(fit2, trend=TRUE)
	markers = marker.fc(fit2, cont.matrix=cont.matrix, log2.threshold = log2(2))
	return(markers)
}

prepare_train<-function(train, train_cell_names){
	# Generate phenodata for reference matrix C
	train_cellID = train
	colnames(train_cellID) = train_cell_names
	# reference matrix (C) + refProfiles.var from TRAINING dataset
	cellType <- colnames(train)
	group = list()
	for(i in unique(cellType)){ 
		group[[i]] <- which(cellType %in% i)
	}
	#C should be made with the mean (not sum) to agree with the way markers were selected
	C = lapply(group,function(x) Matrix::rowMeans(train[,x])) 
	C = round(do.call(cbind.data.frame, C))

	refProfiles.var = lapply(group,function(x) train[,x])
# 	refProfiles.var = lapply(refProfiles.var, function(x) matrixStats::rowSds(Matrix::as.matrix(x)))
    refProfiles.var = lapply(refProfiles.var, function(x) sparseMatrixStats::rowSds(x))   
	refProfiles.var = round(do.call(cbind.data.frame, refProfiles.var))
	rownames(refProfiles.var) <- rownames(train)

	markers = marker_selection(train)
	return(list('train_cellID'=train_cellID, 'C'=C, 'ref' = refProfiles.var, 'markers'=markers))
}

self_reference<-function(param){
	
	if(length(param)!=11){

		print("Please check that all required parameters are indicated or are correct")
		print("Example usage for bulk deconvolution methods: 'Rscript Master_deconvolution.R baron none bulk TMM all nnls 100 none 1'")
		print("Example usage for single-cell deconvolution methods: 'Rscript Master_deconvolution.R baron none sc TMM TMM MuSiC 100 none 1'")
		stop()
	} 

	flag = FALSE

	### arguments
	# 1: dataset name
	# 2: Transformation: none (defalt), log, sqrt, vst
	# 3: deconvolution type: bulk, sc
	# 4: Normalization for C, normalization
	# 5: Normalization for T, marker strategy
	# 6: Deconvolution method
	# 7: number of cells used
	# 8: remove cell type or not (none: default)
	# 9: number of cores used.
	# 10: sampleCT.
	# 11: propsample.

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
	P <- Xtest$P

	#-------------------------------------------------------
	### Transformation, scaling/normalization, marker selection for bulk deconvolution methods and deconvolution:
	if(deconv_type == "bulk"){
#~ 		T = Transformation(Xtest$T, transformation)
#~ 		C = Transformation(Xtrain$C, transformation)

#~ 		T = Scaling(T, normalization)
#~ 		C = Scaling(C, normalization)
		
		T = Scaling(Xtest$T, normalization)
		C = Scaling(Xtrain$C, normalization)
		
		T = Transformation(T, transformation)
		C = Transformation(C, transformation)

		

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

#~ 		T = Transformation(Xtest$T, transformation)
#~ 		C = Transformation(Xtrain$train_cellID, transformation)

#~ 		T = Scaling(T, normalization_scT)
#~ 		C = Scaling(C, normalization_scC)

		T = Scaling(Xtest$T, normalization_scT)
		C = Scaling(Xtrain$train_cellID, normalization_scC)
		
		T = Transformation(T, transformation)
		C = Transformation(C, transformation)

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
	RESULTS = Deconvolution(T = T, C = C, method = method, P = P, elem = to_remove, refProfiles.var = Xtrain$ref, STRING = as.character(sample(1:10000, 1)), marker_distrib = marker_distrib, phenoDataC = pDataC) 
    
#     saveRDS(RESULTS, 'RESULTS.rds')
    
# 	RESULTS = RESULTS %>% dplyr::summarise(RMSE = sqrt(mean((observed_values-expected_values)^2)) %>% round(.,4), 
# 										   Pearson=cor(observed_values,expected_values) %>% round(.,4))
# 	print(RESULTS)
	return(RESULTS)
}

self_reference_pro<-function(param){
	
	if(length(param)!=12){

		print("Please check that all required parameters are indicated or are correct")
		print("Example usage for bulk deconvolution methods: 'Rscript Master_deconvolution.R baron none bulk TMM all nnls 100 none 1'")
		print("Example usage for single-cell deconvolution methods: 'Rscript Master_deconvolution.R baron none sc TMM TMM MuSiC 100 none 1'")
		stop()
	} 

	flag = FALSE

	### arguments
	# 1: dataset name
	# 2: Transformation: none (defalt), log, sqrt, vst
	# 3: deconvolution type: bulk, sc
	# 4: Normalization for C, normalization
	# 5: Normalization for T, marker strategy
	# 6: Deconvolution method
	# 7: number of cells used
	# 8: remove cell type or not (none: default)
	# 9: number of cores used.
	# 10: sampleCT.
	# 11: propsample.
	# 12: Normalize first (T) or Transform first(F). 

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
    NormTrans = param[12]
    
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
	P <- Xtest$P

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
	RESULTS = Deconvolution(T = T, C = C, method = method, P = P, elem = to_remove, refProfiles.var = Xtrain$ref, STRING = as.character(sample(1:10000, 1)), marker_distrib = marker_distrib, phenoDataC = pDataC) 
    
#     saveRDS(RESULTS, 'RESULTS.rds')
    
# 	RESULTS = RESULTS %>% dplyr::summarise(RMSE = sqrt(mean((observed_values-expected_values)^2)) %>% round(.,4), 
# 										   Pearson=cor(observed_values,expected_values) %>% round(.,4))
# 	print(RESULTS)
	return(RESULTS)
}

cross_reference<-function(param){
	
	if(length(param)!=12){

		print("Please check that all required parameters are indicated or are correct")
		print("Example usage for bulk deconvolution methods: 'Rscript Master_deconvolution.R baron none bulk TMM all nnls 100 none 1'")
		print("Example usage for single-cell deconvolution methods: 'Rscript Master_deconvolution.R baron none sc TMM TMM MuSiC 100 none 1'")
		stop()
	} 

	flag = FALSE

	### arguments
	# 1: dataset1 name for reference
	# 2: dataset2 name for pseudobulk
	# 3: Transformation: none (defalt), log, sqrt, vst
	# 4: deconvolution type: bulk, sc
	# 5: Normalization for C, normalization
	# 6: Normalization for T, marker strategy
	# 7: Deconvolution method
	# 8: number of cells used
	# 9: remove cell type or not (none: default)
	# 10: number of cores used.
	# 11: sampleCT.
	# 12: propsample.

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
	P <- Xtest$P

	#-------------------------------------------------------
	### Transformation, scaling/normalization, marker selection for bulk deconvolution methods and deconvolution:
	if(deconv_type == "bulk"){
#~ 		T = Transformation(Xtest$T, transformation)
#~ 		C = Transformation(Xtrain$C, transformation)

#~ 		T = Scaling(T, normalization)
#~ 		C = Scaling(C, normalization)

		T = Scaling(Xtest$T, normalization)
		C = Scaling(Xtrain$C, normalization)
		
		T = Transformation(T, transformation)
		C = Transformation(C, transformation)

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

#~ 		T = Transformation(Xtest$T, transformation)
#~ 		C = Transformation(Xtrain$train_cellID, transformation)

#~ 		T = Scaling(T, normalization_scT)
#~ 		C = Scaling(C, normalization_scC)

		T = Scaling(Xtest$T, normalization_scT)
		C = Scaling(Xtrain$train_cellID, normalization_scC)
		
		T = Transformation(T, transformation)
		C = Transformation(C, transformation)

		#If a cell type is removed, only meaningful mixtures where that CT was present (proportion < 0) are kept:
		if(to_remove != "none"){

			T <- T[,P[to_remove,] != 0]
			C <- C[,pDataC$cellType != to_remove]
			P <- P[!rownames(P) %in% to_remove, colnames(T)]
			pDataC <- pDataC[pDataC$cellType != to_remove,]

		}
		marker_distrib = NULL
	}
	RESULTS = Deconvolution(T = T, C = C, method = method, P = P, elem = to_remove, refProfiles.var = Xtrain$ref, STRING = as.character(sample(1:10000, 1)), marker_distrib = marker_distrib, phenoDataC = pDataC) 
# 	RESULTS = RESULTS %>% dplyr::summarise(RMSE = sqrt(mean((observed_values-expected_values)^2)) %>% round(.,4), 
# 										   Pearson=cor(observed_values,expected_values) %>% round(.,4))
# 	print(RESULTS)
	return(RESULTS)
}

bulk_reference<-function(param){
	
	if(length(param)!=10){

		print("Please check that all required parameters are indicated or are correct")
		print("Example usage for bulk deconvolution methods: 'Rscript Master_deconvolution.R baron none bulk TMM all nnls 100 none 1'")
		print("Example usage for single-cell deconvolution methods: 'Rscript Master_deconvolution.R baron none sc TMM TMM MuSiC 100 none 1'")
		stop()
	} 

	flag = FALSE

	### arguments
	# 1: bulk name
	# 2: reference dataset name
	# 3: Transformation: none (defalt), log, sqrt, vst
	# 4: deconvolution type: bulk, sc
	# 5: Normalization for C, normalization
	# 6: Normalization for T, marker strategy
	# 7: Deconvolution method
	# 8: number of cells used
	# 9: remove cell type or not (none: default)
	# 10: number of cores used.

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

	#-------------------------------------------------------
	### Read single cell data and metadata
	X1 = read_bulk(bulk)
	X2 = read_data(dataset)

	#-------------------------------------------------------
	### QC
	if(FALSE){
#~ 		X1<-QC(X1)
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
#~ 	test <- X2$data
#~ 	colnames(test) <- X2$original_cell_names

#~ 	Xtest <- Generator(sce = test, phenoData = X2$pData, Num.mixtures = 1000, pool.size = number_cells)
#~ 	P <- Xtest$P

	P<-as.data.frame(matrix(0.1, nrow=length(unique(colnames(X2$data))), ncol=dim(X1$data)[2]))
	rownames(P)<-unique(colnames(X2$data))
	colnames(P)<-colnames(X1$data)
	Xtest<-list('T'=X1$data, 'P'=P)
	# Xtest$P
	# Xtest$T

	#-------------------------------------------------------
	### Transformation, scaling/normalization, marker selection for bulk deconvolution methods and deconvolution:
	if(deconv_type == "bulk"){
#~ 		T = Transformation(Xtest$T, transformation)
#~ 		C = Transformation(Xtrain$C, transformation)

#~ 		T = Scaling(T, normalization)
#~ 		C = Scaling(C, normalization)

		T = Scaling(Xtest$T, normalization)
		C = Scaling(Xtrain$C, normalization)
		
		T = Transformation(T, transformation)
		C = Transformation(C, transformation)

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

#~ 		T = Transformation(Xtest$T, transformation)
#~ 		C = Transformation(Xtrain$train_cellID, transformation)

#~ 		T = Scaling(T, normalization_scT)
#~ 		C = Scaling(C, normalization_scC)

		T = Scaling(Xtest$T, normalization_scT)
		C = Scaling(Xtrain$train_cellID, normalization_scC)
		
		T = Transformation(T, transformation)
		C = Transformation(C, transformation)

		#If a cell type is removed, only meaningful mixtures where that CT was present (proportion < 0) are kept:
		if(to_remove != "none"){

			T <- T[,P[to_remove,] != 0]
			C <- C[,pDataC$cellType != to_remove]
			P <- P[!rownames(P) %in% to_remove, colnames(T)]
			pDataC <- pDataC[pDataC$cellType != to_remove,]

		}
		marker_distrib = NULL
	}
	RESULTS = Deconvolution(T = T, C = C, method = method, P = P, elem = to_remove, refProfiles.var = Xtrain$ref, STRING = as.character(sample(1:10000, 1)), marker_distrib = marker_distrib, phenoDataC = pDataC) 
	return(RESULTS)
}

bulk_2references<-function(param){
	
	if(length(param)!=11){

		print("Please check that all required parameters are indicated or are correct")
		print("Example usage for bulk deconvolution methods: 'Rscript Master_deconvolution.R baron none bulk TMM all nnls 100 none 1'")
		print("Example usage for single-cell deconvolution methods: 'Rscript Master_deconvolution.R baron none sc TMM TMM MuSiC 100 none 1'")
		stop()
	} 

	flag = FALSE

	### arguments
	# 1: bulk
	# 2: reference1
	# 3: reference2
	# 4: Transformation: none (defalt), log, sqrt, vst
	# 5: deconvolution type: bulk, sc
	# 6: Normalization for C, normalization
	# 7: Normalization for T, marker strategy
	# 8: Deconvolution method
	# 9: number of cells used
	# 10: remove cell type or not (none: default)
	# 11: number of cores used.
	
	param1<-c(param[1], param[2], param[4], param[5], param[6], param[7], param[8], param[9], param[10], param[11])
	param2<-c(param[1], param[3], param[4], param[5], param[6], param[7], param[8], param[9], param[10], param[11])
	
	RESULTS1 = bulk_reference(param1)
	RESULTS2 = bulk_reference(param2)

# 	saveRDS(RESULTS1, "RESULTS1.rds")
# 	saveRDS(RESULTS2, "RESULTS2.rds")
    
    RESULTS = merge(RESULTS1, RESULTS2, by=c("CT","tissue"))
    
#     RESULTS = RESULTS %>% dplyr::summarise(RMSE = sqrt(mean((observed_values.x-observed_values.y)^2)) %>% round(.,4), 
# 										   Pearson=cor(observed_values.x,observed_values.y) %>% round(.,4))
	return(RESULTS)
}

if(args[1]=='s'){
#~ 	RESULTS = self_reference(args[2:length(args)])
	RESULTS = self_reference_pro(args[2:length(args)])
}else if(args[1]=='c'){
	RESULTS = cross_reference(args[2:length(args)])
}else if(args[1]=='b'){
	RESULTS = bulk_2references(args[2:length(args)])
}
name = getname(args)
saveRDS(RESULTS, paste0("RDS/",name,"rds"))
