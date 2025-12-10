#' Load Baron Human Pancreas Dataset
#'
#' Loads and processes the Baron et al. human pancreas single-cell RNA-seq dataset
#'
#' @return A list containing:
#'   \item{data}{Expression matrix with genes as rows and cells as columns}
#'   \item{celltypes}{Named vector of cell type annotations}
#' @importFrom scRNAseq BaronPancreasData
#' @importFrom SingleCellExperiment counts
#' @examples
#' 
#' # Load Baron human pancreas dataset
#' baron <- hBaronDataset()
#'
#' # Check dimensions
#' dim(baron$data)
#'
#' # View cell type distribution
#' table(baron$celltypes)
#'
#'
#' @export
hBaronDataset <- function() {
  # Load raw data
  baron_h <- scRNAseq::BaronPancreasData(which = "human")

  # Extract and process cell type annotations
  celltype_mappings <- c(
    "activated_stellate" = "activated stellate",
    "quiescent_stellate" = "quiescent stellate",
    "t_cell" = "T cell"
  )

  baron_h_celltype <- baron_h[["label"]]
  names(baron_h_celltype) <- baron_h@colData@rownames

  # Update cell type labels using vectorized operation
  for (old_type in names(celltype_mappings)) {
    baron_h_celltype[baron_h_celltype == old_type] <- celltype_mappings[old_type]
  }

  # Process expression data
  baron_h_data <- as.matrix(SingleCellExperiment::counts(baron_h))
  rownames(baron_h_data) <- gsub(".", "-", rownames(baron_h_data), fixed = TRUE)

  list(
    data = baron_h_data,
    celltypes = baron_h_celltype
  )
}

#' @title Load Baron Human Pancreas Dataset (Deprecated)
#' @description This function has been renamed. Please use \code{hBaronDataset()} instead.
#' @return A list containing the dataset
#' @export
#' @keywords internal
h.baron_dataset <- function() {
  .Deprecated("hBaronDataset")
  hBaronDataset()
}

#' Load Muraro Pancreas Dataset
#'
#' Loads and processes the Muraro et al. pancreas single-cell RNA-seq dataset
#'
#' @return A list containing:
#'   \item{data}{Expression matrix with genes as rows and cells as columns}
#'   \item{celltypes}{Named vector of cell type annotations}
#' @importFrom scRNAseq MuraroPancreasData
#' @examples
#'
#' # Load Muraro pancreas dataset
#' muraro <- muraro_dataset()
#'
#' # Check dataset dimensions
#' dim(muraro$data)
#'
#' # View available cell types
#' table(muraro$celltypes)
#'
#' # Check number of cells per type
#' sort(table(muraro$celltypes), decreasing = TRUE)
#' 
#'
#' @export
muraro_dataset <- function() {
  # Load and filter data
  sceM <- scRNAseq::MuraroPancreasData()
  sceM <- sceM[, !is.na(sceM$label) & sceM$label != "unclear"]

  # Process expression data
  sceM.data <- as.matrix(counts(sceM))
  rownames(sceM.data) <- gsub("\\__chr[0-9]*", "", rownames(sceM.data))
  rownames(sceM.data) <- gsub(".", "-", rownames(sceM.data), fixed = TRUE)

  # Process cell types
  celltype_mappings <- c(
    "pp" = "gamma",
    "duct" = "ductal"
  )

  sceM.celltype <- sceM$label
  for (old_type in names(celltype_mappings)) {
    sceM.celltype[sceM.celltype == old_type] <- celltype_mappings[old_type]
  }

  # Remove mesenchymal cells
  keep_cells <- sceM.celltype != "mesenchymal"
  sceM.data <- sceM.data[, keep_cells]
  sceM.celltype <- sceM.celltype[keep_cells]

  names(sceM.celltype) <- colnames(sceM.data)

  list(
    data = sceM.data,
    celltypes = sceM.celltype
  )
}