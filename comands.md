
# self-reference
Rscript CATD.R ../simdataR/3type_batch1.rds none bulk TMM all nnls 100 none 1

Rscript CATD.R s ../simdataR/100type_batch1.rds none bulk TMM all nnls 100 none 1 F T

Rscript CATD.R c ../simdataR/10type_batch1.rds ../simdataR/10type_batch2.rds none bulk TMM all nnls 10000 none 1 F T

Rscript CATD.R b T.rds example.rds example.rds none bulk TMM all nnls 10000 none 1


git add .
git commit -m "first commit"
git push -u origin main