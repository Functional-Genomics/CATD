
# self-reference
Rscript CATD.R ../simdataR/3type_batch1.rds none bulk TMM all nnls 100 none 1

/nfs/leia/research/ma/chichau/bin/anaconda3/envs/r4-base/bin/Rscript CATD.R s ../simdataR/3type_batch1.rds none bulk TMM all nnls 100 none 1 F T T
/nfs/leia/research/ma/chichau/bin/anaconda3/envs/r4-base/bin/Rscript CATD.R s ../simdataR/3type_batch1.rds none sc TMM all bseqsc 100 none 1 F T T

Rscript CATD.R s ../simdataR/100type_batch1.rds none bulk TMM all nnls 100 none 1 F T

Rscript CATD.R c ../simdataR/10type_batch1.rds ../simdataR/10type_batch2.rds none bulk TMM all nnls 10000 none 1 F T

Rscript CATD.R b T.rds example.rds example.rds none bulk TMM all nnls 10000 none 1

/nfs/leia/research/ma/chichau/bin/anaconda3/envs/r4-base/bin/Rscript CATD_script.R c ../seurat_BRR/Adams2020DS90K.rds ../seurat_BRR/Reyfman2019.rds none bulk column all nnls 10000 none 1 F T T

/nfs/leia/research/ma/chichau/bin/anaconda3/envs/r4-base/bin/Rscript CATD_script.R b ../bulk/Lawlor2017.rds ../seurat_BRR/Lawlor2017.rds ../seurat_BRR/BaronHuman.rds none bulk column all nnls 10000 none 1 T

git add .
git commit -m "Chichau 9Jun2021"
git push -u origin main


echo "# CATD_R_package" >> README.md
git init
git add README.md
git commit -m "first commit"
git branch -M main
git remote add origin https://github.com/Functional-Genomics/CATD_R_package.git
git push -u origin main


R CMD build CATD
Rscript -e 'devtools::install_local("CATD_0.0.1.tar.gz")'


git clone https://github.com/Functional-Genomics/CATD_R_package.git
cd CATD_R_package
conda env create -f CATD.yml
conda activate mycatd
Rscript -e 'devtools::install_github("Functional-Genomics/CATD_R_package", dependencies = TRUE, force = TRUE)'
conda deactivate && conda remove --name mycatd --all

python -m ipykernel install --name myshare --user


sceasy::convertFormat("combined.h5", from="anndata", to="seurat",
                       outFile="combined.rds")

library(Seurat)
ad = readRDS("combined.rds")
ad1 = subset(ad, subset = (SCCAF_Li == "CD16 monocytes"))
Idents(ad1)<-ad1@meta.data[,"CoVID-19 severity"]
de = FindAllMarkers( ad1, test.use = 'MAST')

de %>%
    group_by(cluster) %>%
    top_n(n = 2, wt = avg_log2FC)


conda create --name envr   ###create the environment "envr"
conda activate envr       ### activate the environment ##then you install the following things 
conda install -c conda-forge tensorflow
conda install -c anaconda ipykernel
python -m ipykernel install --user --name=myscanpy
###then start the jupyter notebook you will be able to see the conda environments



1. Create a skeleton package
usethis::create_package(
"pkgname"
)

2. License code
usethis::use_mit_license()

3. Enable code documentation inline
usethis::use_roxygen_md()

4. Create an internal package-level help documentation file
usethis::use_package_doc()

5. Add R code file our package's  R/ folder 

usethis::use_r(
"test"
) 


devtools::document()

#' @title
#' 
#' @description
#' 
#' @details
#' 
#' @param 
#' @param 
#' @param 
#' @param 
#' @param 
#' 
#' @return
#' 
#' @example
#' 


CATD.R: framework.R + CATD_script.R
helper_functions.R: benchmark.R
new: result.R + format.R



Rscript CATD_script.R s data/3type_batch1.rds none sc TMM all BisqueRNA 100 none 1 F T T
Rscript CATD_script.R s data/3type_batch1.rds none sc TMM all CPM 100 none 1 F T T
Rscript CATD_script.R s data/3type_batch1.rds none sc TMM all MuSiC 100 none 1 F T T 
Rscript CATD_script.R s data/3type_batch1.rds none sc TMM all SCDC 100 none 1 F T T 
Rscript CATD_script.R s data/3type_batch1.rds none sc TMM all bseqsc 100 none 1 F T T 
Rscript CATD_script.R s data/3type_batch1.rds none sc TMM all deconvSeq 100 none 1 F T T
Rscript CATD_script.R s data/3type_batch1.rds none sc TMM all DWLS 100 none 1 F T T
Rscript CATD_script.R s data/3type_batch1.rds none sc TMM all TIMER 100 none 1 F T T 

Rscript CATD_script.R s data/3type_batch1.rds none bulk TMM all CAMmarker 100 none 1 F T T
Rscript CATD_script.R s data/3type_batch1.rds none bulk TMM all CDSeq 100 none 1 F T T 
Rscript CATD_script.R s data/3type_batch1.rds none bulk TMM all CIBERSORT 100 none 1 F T T
Rscript CATD_script.R s data/3type_batch1.rds none bulk TMM all DCQ 100 none 1 F T T
Rscript CATD_script.R s data/3type_batch1.rds none bulk TMM all DSA 100 none 1 F T T
Rscript CATD_script.R s data/3type_batch1.rds none bulk TMM all DeconRNASeq 100 none 1 F T T
Rscript CATD_script.R s data/3type_batch1.rds none bulk TMM all EPIC 100 none 1 F T T
Rscript CATD_script.R s data/3type_batch1.rds none bulk TMM all EpiDISH 100 none 1 F T T
Rscript CATD_script.R s data/3type_batch1.rds none bulk TMM all FARDEEP 100 none 1 F T T
Rscript CATD_script.R s data/3type_batch1.rds none bulk TMM all OLS 100 none 1 F T T
Rscript CATD_script.R s data/3type_batch1.rds none bulk TMM all RLR 100 none 1 F T T
Rscript CATD_script.R s data/3type_batch1.rds none bulk TMM all deconf 100 none 1 F T T
Rscript CATD_script.R s data/3type_batch1.rds none bulk TMM all dtangle 100 none 1 F T T
Rscript CATD_script.R s data/3type_batch1.rds none bulk TMM all elasticNet 100 none 1 F T T
Rscript CATD_script.R s data/3type_batch1.rds none bulk TMM all lasso 100 none 1 F T T
Rscript CATD_script.R s data/3type_batch1.rds none bulk TMM all ridge 100 none 1 F T T
Rscript CATD_script.R s data/3type_batch1.rds none bulk TMM all nnls 100 none 1 F T T
Rscript CATD_script.R s data/3type_batch1.rds none bulk TMM all proportionsInAdmixture 100 none 1 F T T
Rscript CATD_script.R s data/3type_batch1.rds none bulk TMM all ssFrobenius 100 none 1 F T T
Rscript CATD_script.R s data/3type_batch1.rds none bulk TMM all ssKL 100 none 1 F T T



deconf??


git add .
git commit -m "x"
git push -u origin main
