#' Perform Sub-clustering Analysis on CSFNMF Results
#'
#' @description
#' Identifies and predicts cell sub-clusters based on the CSFNMF factorization results.
#' This function performs sub-clustering on both reference and query datasets,
#' preserving the hierarchical relationship between main cell types and sub-clusters.
#'
#' @param object A CSFNMF object with completed training
#' @param resolution Clustering resolution parameter (default: 0.5)
#' @param algorithm Clustering algorithm (1: Louvain, 2: Leiden, 3: SLM, 4: LabelProp)
#' @param min_cells Minimum number of cells required for sub-clustering (default: 30)
#' @param verbose Show progress messages (default: TRUE)
#'
#' @return Updated CSFNMF object with sub-clustering results stored in annotations
#'
#' @details
#' The function performs the following steps:
#' 1. Identifies sub-clusters in the reference dataset using Seurat clustering
#' 2. Updates reference annotations with sub-cluster information
#' 3. Predicts sub-clusters for the query dataset using SingleR
#' 4. Stores both original and sub-clustered annotations
#'
#' @examples
#' \dontrun{
#' # Train CSFNMF model
#' object_trained <- RunCSFNMF(object)
#'
#' # Perform sub-clustering with default parameters
#' object_subclustered <- SubClusCSFNMF(object_trained)
#'
#' # Use different resolution and algorithm
#' object_subclustered <- SubClusCSFNMF(
#'   object_trained,
#'   resolution = 0.8,
#'   algorithm = 2  # Use Leiden algorithm
#' )
#'
#' # Access sub-clustering results
#' head(object_subclustered@annotation$subclusters)
#' }
#'
#' @export
SubClusCSFNMF <- function(object,
                          resolution = 0.5,
                          algorithm = 1,
                          min_cells = 30,
                          verbose = TRUE) {
  
  # Input validation
  if (!inherits(object, "csfnmf")) {
    stop("Input must be a CSFNMF object")
  }
  if (!resolution > 0) {
    stop("Resolution must be positive")
  }
  if (!algorithm %in% 1:4) {
    stop("Algorithm must be between 1 and 4")
  }
  if (!min_cells > 0) {
    stop("min_cells must be positive")
  }
  
  report <- create_reporter(verbose)
  
  # Find sub-clusters in reference data
  report("Finding sub-clusters in reference data")
  ref_subclusters <- find_subclusters(
    object@train_object,
    resolution = resolution,
    algorithm = algorithm,
    min_cells = min_cells,
    verbose = verbose
  )
  
  # Update reference annotations
  report("Updating reference annotations")
  object@train_object@annotation["original_celltype"] <- 
    object@train_object@annotation$celltype
  object@train_object@annotation$celltype <- 
    ref_subclusters[rownames(object@train_object@annotation)]
  
  # Reorder and update object
  object@train_object <- reorder_data(object@train_object, "celltype")
  object@train_object@annotation["celltype.code"] <- 
    encode_celltypes(object@train_object@annotation$celltype)
  object@train_object@H <- 
    object@train_object@H[, rownames(object@train_object@annotation)]
  object@train_object <- calculate_help_matrices(object@train_object)
  
  # Predict sub-clusters for validation data
  report("Predicting sub-clusters for query data")
  valid_subclusters <- predict_subclusters(
    object,
    verbose = verbose
  )
  
  # Update validation annotations
  object@train_object@test.annotation["original_celltype"] <- 
    object@train_object@test.annotation$celltype
  object@train_object@test.annotation$celltype <- 
    valid_subclusters[rownames(object@train_object@test.annotation)]
  
  # Store all sub-clusters
  object@annotation["subclusters"] <- 
    c(ref_subclusters, valid_subclusters)[rownames(object@annotation)]
  
  report("Sub-clustering complete")
  object
}

#' Find Sub-clusters in Reference Data
#'
#' @param train_object Training object
#' @param resolution Clustering resolution
#' @param algorithm Clustering algorithm
#' @param min_cells Minimum cells per cluster
#' @param verbose Show progress
#' @return Vector of sub-cluster assignments
#' @importFrom Seurat CreateSeuratObject FindVariableFeatures FindNeighbors FindClusters
#' @keywords internal
find_subclusters <- function(train_object,
                             resolution,
                             algorithm,
                             min_cells,
                             verbose = TRUE) {
  # Get cell types
  celltype <- train_object@annotation$celltype
  names(celltype) <- rownames(train_object@annotation)
  new_celltype <- celltype
  
  # Process each cell type
  for (type in unique(celltype)) {
    cells_of_type <- which(celltype == type)
    n_cells <- length(cells_of_type)
    
    if (n_cells > min_cells) {
      # Extract data for current type
      data_subset <- train_object@H[, cells_of_type, drop = FALSE]
      
      # Create and process Seurat object
      seurat_obj <- tryCatch({
        seurat_object <- CreateSeuratObject(counts = data_subset)
        assay_v3 <- CreateAssayObject(
          counts = seurat_object[["RNA"]]$counts
        )
        
        seurat_object[["RNA"]] <- assay_v3
        seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
        seurat_object <- FindNeighbors(seurat_object, dims = NULL, annoy.metric="cosine", verbose = FALSE)
        seurat_object <- FindClusters(seurat_object, resolution = resolution, algorithm = algorithm, verbose = FALSE)
        seurat_object
      }, error = function(e) {
        warning(sprintf("Failed to cluster %s: %s", type, e$message))
        NULL
      })
      
      # Update cell type labels if clustering successful
      if (!is.null(seurat_obj)) {
        clusters <- Seurat::Idents(seurat_obj)
        if (length(unique(clusters)) > 1) {
          new_celltype[names(clusters)] <- paste0(type, "_",
                                                  as.numeric(clusters))
        }
      }
    }
  }
  
  new_celltype
}

#' Predict Sub-clusters for Query Data
#'
#' @param object CSFNMF object
#' @param verbose Show progress
#' @return Vector of predicted sub-clusters
#' @importFrom SingleR SingleR
#' @keywords internal
predict_subclusters <- function(object, verbose = TRUE) {
  # Initialize
  celltype <- object@train_object@test.annotation$celltype
  names(celltype) <- rownames(object@train_object@test.annotation)
  new_subcluster <- celltype
  
  # Get reference data
  ref_subcluster <- object@train_object@annotation$celltype
  names(ref_subcluster) <- rownames(object@train_object@annotation)
  
  h_ref <- object@H
  h_project <- object@H[, names(new_subcluster)]
  
  # Process each cell type
  for (type in unique(celltype)) {
    # Get cells of current type
    type_cells <- names(celltype)[celltype == type]
    type_h_project <- as.data.frame(h_project[, type_cells, drop = FALSE])
    
    # Get reference sub-clusters
    ref_cells <- which(gsub('_[0-9]*', '', ref_subcluster) == type)
    type_ref_subcluster <- ref_subcluster[ref_cells]
    
    if (length(unique(type_ref_subcluster)) > 1) {
      # Predict sub-clusters
      type_h_ref <- as.data.frame(h_ref[, names(type_ref_subcluster),
                                        drop = FALSE])
      
      pred <- SingleR(
        test = type_h_project,
        ref = type_h_ref,
        labels = type_ref_subcluster,
        de.n = 50
      )
      
      new_subcluster[pred@rownames] <- pred@listData[["labels"]]
    }
  }
  
  new_subcluster
}