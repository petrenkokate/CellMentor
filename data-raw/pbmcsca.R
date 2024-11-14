# data-raw/pbmcsca.R
# library(Seurat)
# library(SeuratData)
# InstallData("pbmcsca")
# data("pbmcsca")
# pbmcsca <- UpdateSeuratObject(pbmcsca)
# 
# saveRDS(pbmcsca, file='~/CellMentor/data-raw/pbmcsca_raw.rds')
# print('DONE!')
pbmcsca <- readRDS('~/CellMentor/data-raw/pbmcsca_raw.rds')
usethis::use_data(pbmcsca, overwrite = TRUE)
