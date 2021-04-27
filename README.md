# CATD: Critical Assessment of Transcriptomics Deconvolution

This package extends the functionalities in [deconv_benchmark](https://github.com/favilaco/deconv_benchmark), which was published at: [Nature Communications](https://doi.org/10.1038/s41467-020-19015-1). 

# Citation 
Zhichao Miao, Anna Vathrakokili Pournara, Irene Papatheodoru. Multiple-level benchmarking of cell type deconvolution. 


# Installation
The scripts in CATD include Master_deconvolution.R and helper_functions.R, which do not require an installation. 
The installation instruction for the deconvolution programs can be found at [install.md](https://github.com/Functional-Genomics/CATD/blob/main/install.md). 

# Input format
The input data are formatted as [ExpressionSet](https://www.bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf).
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


# Experimental design


# Example command

# 

