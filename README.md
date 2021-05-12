<pre>
<img src="https://github.com/Functional-Genomics/CATD/blob/main/CATD.jpg" alt="drawing" width="100"/>

  /$$$$$$   /$$$$$$  /$$$$$$$$ /$$$$$$$ 
 /$$__  $$ /$$__  $$|__  $$__/| $$__  $$
| $$  \__/| $$  \ $$   | $$   | $$  \ $$
| $$      | $$$$$$$$   | $$   | $$  | $$
| $$      | $$__  $$   | $$   | $$  | $$
| $$    $$| $$  | $$   | $$   | $$  | $$
|  $$$$$$/| $$  | $$   | $$   | $$$$$$$/
 \______/ |__/  |__/   |__/   |_______/ 

</pre>                                 


# CATD: Critical Assessment of Transcriptomics Deconvolution

This package extends the functionalities in [deconv_benchmark](https://github.com/favilaco/deconv_benchmark), which was published at: [Nature Communications](https://doi.org/10.1038/s41467-020-19015-1). 

# Citation 
Zhichao Miao, Anna Vathrakokili Pournara, Irene Papatheodoru. Multiple-level benchmarking of cell type deconvolution. 


# Installation
The scripts in CATD include CATD.R and helper_functions.R, which do not require an installation. 
The installation instruction for the deconvolution programs can be found at [install.md](https://github.com/Functional-Genomics/CATD/blob/main/install.md). 

# Input format
The input reference data (a single-cell dataset) is formatted as [ExpressionSet](https://www.bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf).
```
ExpressionSet (storageMode: lockedEnvironment)
assayData: 80 features, 1000 samples 
  element names: exprs 
protocolData: none
phenoData
  sampleNames: cell_1 cell_2 ... cell_1000 (1000 total)
  varLabels: cellID cellType sampleID
  varMetadata: labelDescription
featureData: none
experimentData: use 'experimentData(object)'
Annotation:  
```
The expression matrix is available at: `x@assayData$exprs`. Rows for genes, columns for cells. `A matrix: M × N of type int`. 
The metadata is available at: `x@phenoData@data`, which includes THREE columns: cellID, cellType and sampleID.

An example below:
```
A data.frame: 1000 × 3
		cellID	cellType	sampleID
		<chr>	<chr>	<chr>
cell_1	cell_1	cell_type_1	indiv1
cell_2	cell_2	cell_type_1	indiv2
cell_3	cell_3	cell_type_1	indiv1
cell_4	cell_4	cell_type_1	indiv2
cell_5	cell_5	cell_type_1	indiv1
cell_6	cell_6	cell_type_1	indiv2
cell_7	cell_7	cell_type_1	indiv1
cell_8	cell_8	cell_type_1	indiv2
cell_9	cell_9	cell_type_1	indiv1
cell_10	cell_10	cell_type_1	indiv2
......
```

## Format transformation from python
Python formatted [AnnData](https://anndata.readthedocs.io/en/latest/anndata.AnnData.html) can be transformed to required ExpressionSet format using the `anndata2expressionset` function in [helper_functions.R](https://github.com/Functional-Genomics/CATD/blob/main/helper_functions.R). 

An example here: `anndata2expressionset('downsample200/Ximerakis2019.h5', 'downsample200R/Ximerakis2019.rds')`. 


# Experimental design
In this evaluation framework, we have THREE levels of experiments: 1. self-reference, 2. cross-reference and 3. real bulk. 
## 1. self-reference
Both the pseudobulk and the cell type reference are generated from the same dataset, which is the same setting as the [Nature Communications](https://doi.org/10.1038/s41467-020-19015-1) paper. 

In this case, the `self_reference` function in [CATD.R](https://github.com/Functional-Genomics/CATD/blob/main/CATD.R) is used. 

## 2. cross-refernce
When we have two reference (single-cell) datasets for the same type of tissue available, one of the dataset can be used to generate pseudobulk and the other one can be used as cell type reference. 

In this case, the `cross_reference` function in [CATD.R](https://github.com/Functional-Genomics/CATD/blob/main/CATD.R) is used. 

## 3. real bulk
When we have real bulk RNA-seq data available, we can deconvolve real bulk data with a single-cell reference dataset. However, there is no way to know the cell type proportion in the bulk dataset in priori. Therefore, we use TWO reference datasets to deconvolve the real bulk dataset. 
We assume that a program can deconvolve the real bulk dataset with the reference dataset very well (suppose the accuracy is 100%). The resulted cell proportion can be used as a gold standard to compare with the deconvolution using the second reference dataset. There can be possibility that the deconvolution with the two reference datasets are accidentally equally wrong. However, such a situation can be almost excluded after applying the self-reference and cross-reference evaluations. 

The `bulk_reference` function in [CATD.R](https://github.com/Functional-Genomics/CATD/blob/main/CATD.R) can deconvolve a real bulk dataset with a reference dataset, while the `bulk_2references` function compares the deconvolutions of a real bulk dataset using two reference datasets. 

# Example command

> 1. self-reference
```
# With the example we provided with this repository + no cell type removed:
Rscript CATD.R s example.rds none bulk TMM all nnls 100 none 1 F T
	#Expected output:
	#        RMSE   Pearson
	#1     0.0351    0.9866
```

> 2. cross-reference
```
# With the example we provided with this repository + no cell type removed:
Rscript CATD.R c example1.rds example2.rds none bulk TMM all nnls 10000 none 1 F T
	#Expected output:
	#        RMSE   Pearson
	#1     0.0351    0.9866
```

> 3. real bulk
```
# With the example we provided with this repository + no cell type removed:
Rscript CATD.R b bulk.rds example1.rds example2.rds none bulk TMM all nnls 10000 none 1
	#Expected output:
	#        RMSE   Pearson
	#1     0.0351    0.9866
```


