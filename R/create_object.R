#' Create CSFNMF Object for Cell Type Analysis
#'
#' @description
#' Creates and initializes a Constrained Supervised Factorization NMF (CSFNMF) object
#' for analyzing single-cell RNA sequencing data. This is the main function for
#' starting analysis with CellMentor.
#'
#' @param ref_matrix Reference matrix (genes × cells) with known cell types
#' @param ref_celltype Vector of cell type labels for reference cells
#' @param data_matrix Query matrix (genes × cells) to be analyzed
#' @param norm Logical: perform normalization (default: TRUE)
#' @param most.variable Logical: select variable genes (default: TRUE)
#' @param scale Logical: perform scaling (default: TRUE)
#' @param scale_by Character: scaling method, either "cells" or "genes" (default: "cells")
#' @param gene_list Optional vector of genes to include (default: NULL)
#' @param verbose Logical: show progress messages (default: TRUE)
#' @param num_cores Integer: number of cores for parallel processing (default: 1)
#'
#' @return A CSFNMF object containing processed data and annotations
#' 
#' @importFrom methods new
#' @importFrom Matrix rowSums colSums
#' @importFrom SingleR trainSingleR
#' @importFrom Seurat LogNormalize
#' 
#' @examples
#' # Create CSFNMF object with default parameters
#' object <- CreateCSFNMFobject(
#'   ref_matrix = reference_data,
#'   ref_celltype = cell_types,
#'   data_matrix = query_data
#' )
#' 
#' # Create object with custom parameters
#' object <- CreateCSFNMFobject(
#'   ref_matrix = reference_data,
#'   ref_celltype = cell_types,
#'   data_matrix = query_data,
#'   scale_by = "genes",
#'   num_cores = 4
#' )
#' @export
CreateCSFNMFobject <- function(ref_matrix,
                               ref_celltype,
                               data_matrix,
                               norm = TRUE,
                               most.variable = TRUE,
                               scale = TRUE,
                               scale_by = "cells",
                               gene_list = NULL,
                               verbose = TRUE,
                               num_cores = 1) {
  
  # Initialize progress reporting
  report <- create_reporter(verbose)
  report("Starting CSFNMF object creation")
  
  # Monitor memory usage
  with_memory_check({
    # Validate inputs
    report("Validating inputs")
    params <- list(
      scale_by = scale_by,
      num_cores = num_cores
    )
    validate_inputs(ref_matrix, ref_celltype, data_matrix, params)
    
    # Create initial object
    report("Creating CSFNMF object")
    object <- tryCatch({
      methods::new("csfnmf")
    }, error = function(e) {
      stop("Failed to create CSFNMF object: ", e$message)
    })
    
    # Convert matrices to sparse format
    report("Converting matrices to sparse format")
    object@count.matrices <- object@origin.matrices <- tryCatch({
      methods::new(
        "RefDataList",
        ref = to_sparse(ref_matrix),
        data = to_sparse(data_matrix)
      )
    }, error = function(e) {
      stop("Failed to convert matrices to sparse format: ", e$message)
    })
    
    # Set up annotation
    report("Setting up annotations")
    tryCatch({
      names(ref_celltype) <- colnames(object@count.matrices@ref)
      object@annotation <- data.frame(
        celltype = ref_celltype,
        row.names = names(ref_celltype),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      stop("Failed to set up annotations: ", e$message)
    })
    
    # Clean matrices
    report("Cleaning matrices")
    tryCatch({
      # Clean reference matrix
      ref_cleaned <- clean_matrix_rows(object@count.matrices@ref)
      data_cleaned <- clean_matrix_rows(object@count.matrices@data)
      
      if (length(ref_cleaned$removed_rows) > 0) {
        report(sprintf("Removed %d empty rows from reference matrix", 
                       length(ref_cleaned$removed_rows)))
      }
      if (length(data_cleaned$removed_rows) > 0) {
        report(sprintf("Removed %d empty rows from query matrix", 
                       length(data_cleaned$removed_rows)))
      }
      
      object@count.matrices@ref <- ref_cleaned$matrix
      object@count.matrices@data <- data_cleaned$matrix
      
    }, error = function(e) {
      stop("Failed to clean matrices: ", e$message)
    })
    
    # Find common genes
    report("Finding common genes")
    tryCatch({
      common_genes <- find_common_genes(
        object@count.matrices@ref,
        object@count.matrices@data
      )
      
      report(sprintf("Found %d common genes", length(common_genes$genes)))
      
      object@count.matrices@ref <- object@count.matrices@ref[common_genes$genes, ]
      object@count.matrices@data <- object@count.matrices@data[common_genes$genes, ]
      
    }, error = function(e) {
      stop("Failed to find common genes: ", e$message)
    })
    
    # Normalize if requested
    if (norm) {
      report("Normalizing data")
      tryCatch({
        norm_matrices <- list(
          ref = normalize_matrix(object@count.matrices@ref),
          data = normalize_matrix(object@count.matrices@data)
        )
        
        object@count.matrices <- object@norm.matrices <- methods::new(
          "RefDataList",
          ref = norm_matrices$ref,
          data = norm_matrices$data
        )
        
      }, error = function(e) {
        stop("Failed to normalize data: ", e$message)
      })
    }
    
    # Select variable genes if requested
    if (most.variable) {
      report("Selecting variable genes")
      tryCatch({
        var_matrices <- select_variable_genes(object, gene_list, num_cores)
        report(sprintf("Selected %d variable genes", 
                       nrow(var_matrices$ref)))
        
        object@count.matrices <- object@f.select.matrices <- methods::new(
          "RefDataList",
          ref = var_matrices$ref,
          data = var_matrices$data
        )
        
      }, error = function(e) {
        stop("Failed to select variable genes: ", e$message)
      })
    }
    
    # Scale if requested
    if (scale) {
      report(sprintf("Scaling data by %s", scale_by))
      tryCatch({
        scale_fn <- if(scale_by == "cells") {
          scale_by_cells
        } else {
          scale_by_genes
        }
        
        object@count.matrices <- object@scale.matrices <- methods::new(
          "RefDataList",
          ref = scale_fn(object@count.matrices@ref),
          data = scale_fn(object@count.matrices@data)
        )
        
      }, error = function(e) {
        stop("Failed to scale data: ", e$message)
      })
    }
    
    # Encode cell types and reorder
    report("Encoding cell types")
    tryCatch({
      object@annotation$celltype.code <- encode_celltypes(object@annotation$celltype)
      
      report("Reordering data")
      object <- reorder_data(object, "celltype")
      
    }, error = function(e) {
      stop("Failed to encode cell types and reorder data: ", e$message)
    })
    
    # Final validation
    report("Validating final object")
    if (!validate_final_object(object)) {
      stop("Final object validation failed")
    }
    
    report("CSFNMF object creation complete")
    return(object)
    
  }, threshold = 2000)  # Memory check threshold
}

#' Validate Final CSFNMF Object
#'
#' @param object CSFNMF object to validate
#' @return Logical indicating if object is valid
#' @keywords internal
validate_final_object <- function(object) {
  # Check matrix dimensions
  if (ncol(object@count.matrices@ref) != nrow(object@annotation)) {
    return(FALSE)
  }
  
  # Check gene names consistency
  if (!identical(rownames(object@count.matrices@ref),
                 rownames(object@count.matrices@data))) {
    return(FALSE)
  }
  
  # Check annotation completeness
  if (!all(c("celltype", "celltype.code") %in% names(object@annotation))) {
    return(FALSE)
  }
  
  return(TRUE)
}