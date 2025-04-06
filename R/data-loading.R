#' Load PBMC Dataset from SeuratData
#' 
#' @return A list containing processed data and cell types
#' @importFrom SeuratData InstallData
#' @importFrom Seurat UpdateSeuratObject PercentageFeatureSet
#' @export
load_pbmcsca_data <- function() {
  if (!requireNamespace("SeuratData", quietly = TRUE)) {
    message("Installing SeuratData...")
    remotes::install_github("satijalab/seurat-data")
  }
  
  if (!requireNamespace("pbmcsca", quietly = TRUE)) {
    SeuratData::InstallData("pbmcsca")
  }
  library(pbmcsca.SeuratData)
  data("pbmcsca")
  pbmcsca <- UpdateSeuratObject(pbmcsca)
  
  # Process data as before
  pbmcsca[["percent.mt"]] <- PercentageFeatureSet(pbmcsca, pattern = "^MT-")
  pbmcsca <- subset(pbmcsca,
                    subset = nFeature_RNA < 4000 & 
                      percent.mt < 5)
  
  # Extract data and cell types
  pbmcsca.data <- pbmcsca[["RNA"]]@counts
  pbmcsca.celltype <- pbmcsca$CellType
  names(pbmcsca.celltype) <- colnames(pbmcsca.data)
  
  list(
    data = pbmcsca.data,
    celltypes = pbmcsca.celltype
  )
}