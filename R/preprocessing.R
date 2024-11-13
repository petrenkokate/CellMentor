#' @importFrom Matrix rowSums colSums
#' @importFrom SingleR trainSingleR
#' @importFrom methods as
NULL

#' Clean Matrix Rows
#'
#' @description
#' Removes empty rows and replaces NA values with zeros efficiently
#'
#' @param matrix Input sparse matrix
#' @return List with cleaned matrix and removed indices
#' @keywords internal
clean_matrix_rows <- function(matrix) {
  # Replace NAs with zeros using direct assignment
  matrix@x[is.na(matrix@x)] <- 0
  
  # Find zero rows using rowSums
  zero_rows <- which(Matrix::rowSums(matrix) == 0)
  
  # Remove zero rows if any exist
  if (length(zero_rows)) {
    matrix <- matrix[-zero_rows, , drop = FALSE]
  }
  
  list(
    matrix = matrix,
    removed_rows = zero_rows
  )
}

#' Find Common Genes
#'
#' @description
#' Identifies common genes between reference and query datasets
#'
#' @param ref_matrix Reference expression matrix
#' @param data_matrix Query expression matrix
#' @return List of indices and shared genes
#' @keywords internal
find_common_genes <- function(ref_matrix, data_matrix) {
  # Get gene names
  ref_genes <- rownames(ref_matrix)
  data_genes <- rownames(data_matrix)
  
  # Find shared genes
  shared_genes <- intersect(ref_genes, data_genes)
  
  # Get indices using match
  list(
    ref.index = match(shared_genes, ref_genes),
    data.index = match(shared_genes, data_genes),
    genes = shared_genes
  )
}

#' Select Variable Genes
#'
#' @description
#' Identifies highly variable genes using SingleR
#'
#' @param object CSFNMF object
#' @param gene_list Additional genes to include
#' @param num_cores Number of cores for parallel processing
#' @return List of updated matrices
#' @keywords internal
select_variable_genes <- function(object, gene_list = NULL, num_cores = 1) {
  # Get reference data and annotations
  ref_matrix <- object@count.matrices@ref
  annotations <- object@annotation$celltype
  
  # Run SingleR (parallel if multiple cores)
  var_genes_result <- if (num_cores > 1) {
    num_cores <- min(num_cores, parallel::detectCores())
    cl <- parallel::makeCluster(num_cores)
    on.exit(parallel::stopCluster(cl))
    parallel::parLapply(cl, 1L, function(x) {
      SingleR::trainSingleR(ref_matrix, annotations)
    })[[1]]
  } else {
    SingleR::trainSingleR(ref_matrix, annotations)
  }
  
  # Extract variable genes
  variable_genes <- unique(var_genes_result$markers$unique)
  
  # Add additional genes if provided
  if (length(gene_list)) {
    valid_additional <- gene_list[gene_list %in% rownames(ref_matrix)]
    variable_genes <- unique(c(variable_genes, valid_additional))
  }
  
  # Subset matrices
  list(
    ref = ref_matrix[variable_genes, , drop = FALSE],
    data = object@count.matrices@data[variable_genes, , drop = FALSE]
  )
}

#' Encode Cell Types
#'
#' @description
#' Converts cell type labels to numeric codes efficiently
#'
#' @param celltype_list Vector of cell type labels
#' @return Factor with numeric levels
#' @keywords internal
encode_celltypes <- function(celltype_list) {
  # Convert to factor and get numeric codes
  f <- factor(celltype_list)
  codes <- as.numeric(f)
  
  # Store mapping in attributes
  attr(codes, "levels") <- levels(f)
  attr(codes, "mapping") <- setNames(seq_along(levels(f)), levels(f))
  
  codes
}

#' Reorder Reference Data
#'
#' @description
#' Reorders reference data by cell types efficiently
#'
#' @param object CSFNMF object
#' @param reorder_by Column to order by
#' @return Updated object
#' @keywords internal
reorder_data <- function(object, reorder_by) {
  # Get order
  new_order <- order(object@annotation[[reorder_by]])
  
  # Reorder both annotation and matrix columns
  object@annotation <- object@annotation[new_order, , drop = FALSE]
  object@count.matrices@ref <- object@count.matrices@ref[, new_order, drop = FALSE]
  
  object
}