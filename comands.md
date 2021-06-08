
# self-reference
Rscript CATD.R ../simdataR/3type_batch1.rds none bulk TMM all nnls 100 none 1

/nfs/leia/research/ma/chichau/bin/anaconda3/envs/r4-base/bin/Rscript CATD.R s ../simdataR/3type_batch1.rds none bulk TMM all nnls 100 none 1 F T T
/nfs/leia/research/ma/chichau/bin/anaconda3/envs/r4-base/bin/Rscript CATD.R s ../simdataR/3type_batch1.rds none sc TMM all bseqsc 100 none 1 F T T

Rscript CATD.R s ../simdataR/100type_batch1.rds none bulk TMM all nnls 100 none 1 F T

Rscript CATD.R c ../simdataR/10type_batch1.rds ../simdataR/10type_batch2.rds none bulk TMM all nnls 10000 none 1 F T

Rscript CATD.R b T.rds example.rds example.rds none bulk TMM all nnls 10000 none 1


git add .
git commit -m "Chichau 8Jun2021"
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




devtools::document()
