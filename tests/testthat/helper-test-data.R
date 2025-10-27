library(Matrix)

#' Create small test matrices for fast testing
#' @return List with ref_matrix, ref_celltypes, query_matrix, query_celltypes
create_small_test_data <- function(n_genes = 100, 
                                   n_ref_cells = 50, 
                                   n_query_cells = 30,
                                   n_celltypes = 3) {
  
  set.seed(123)  # For reproducibility
  
  # Generate synthetic count data with some structure
  # Cell type 1: high expression in genes 1-30
  # Cell type 2: high expression in genes 31-60
  # Cell type 3: high expression in genes 61-90
  
  ref_matrix <- Matrix::Matrix(0, nrow = n_genes, ncol = n_ref_cells, sparse = TRUE)
  query_matrix <- Matrix::Matrix(0, nrow = n_genes, ncol = n_query_cells, sparse = TRUE)
  
  # Define cell types
  cells_per_type_ref <- n_ref_cells %/% n_celltypes
  cells_per_type_query <- n_query_cells %/% n_celltypes
  
  ref_celltypes <- rep(paste0("Type", 1:n_celltypes), 
                       c(rep(cells_per_type_ref, n_celltypes - 1), 
                         n_ref_cells - (cells_per_type_ref * (n_celltypes - 1))))
  names(ref_celltypes) <- paste0("RefCell_", seq_len(n_ref_cells))
  
  query_celltypes <- rep(paste0("Type", 1:n_celltypes), 
                         c(rep(cells_per_type_query, n_celltypes - 1), 
                           n_query_cells - (cells_per_type_query * (n_celltypes - 1))))
  names(query_celltypes) <- paste0("QueryCell_", seq_len(n_query_cells))
  
  # Generate data with cell-type specific patterns
  genes_per_type <- n_genes %/% n_celltypes
  
  for (i in seq_len(n_celltypes)) {
    # Gene indices for this cell type
    gene_start <- (i - 1) * genes_per_type + 1
    gene_end <- min(i * genes_per_type, n_genes)
    
    # Cell indices for this cell type
    ref_cells_idx <- which(ref_celltypes == paste0("Type", i))
    query_cells_idx <- which(query_celltypes == paste0("Type", i))
    
    # Add high expression for type-specific genes
    ref_matrix[gene_start:gene_end, ref_cells_idx] <- 
      rpois(length(gene_start:gene_end) * length(ref_cells_idx), lambda = 50)
    
    query_matrix[gene_start:gene_end, query_cells_idx] <- 
      rpois(length(gene_start:gene_end) * length(query_cells_idx), lambda = 50)
    
    # Add low background expression for other genes
    ref_matrix[-c(gene_start:gene_end), ref_cells_idx] <- 
      rpois(length(-c(gene_start:gene_end)) * length(ref_cells_idx), lambda = 5)
    
    query_matrix[-c(gene_start:gene_end), query_cells_idx] <- 
      rpois(length(-c(gene_start:gene_end)) * length(query_cells_idx), lambda = 5)
  }
  
  # Add row and column names
  rownames(ref_matrix) <- paste0("Gene_", seq_len(n_genes))
  colnames(ref_matrix) <- names(ref_celltypes)
  
  rownames(query_matrix) <- paste0("Gene_", seq_len(n_genes))
  colnames(query_matrix) <- names(query_celltypes)
  
  list(
    ref_matrix = ref_matrix,
    ref_celltypes = ref_celltypes,
    query_matrix = query_matrix,
    query_celltypes = query_celltypes
  )
}

#' Create minimal test data (for very fast tests)
create_minimal_test_data <- function() {
  create_small_test_data(
    n_genes = 50,
    n_ref_cells = 20,
    n_query_cells = 10,
    n_celltypes = 2
  )
}

#' Create test data with edge cases
create_edge_case_data <- function(case = "single_celltype") {
  
  set.seed(456)
  
  switch(case,
         "single_celltype" = {
           # Only one cell type
           data <- create_small_test_data(n_celltypes = 1)
           data
         },
         "unbalanced" = {
           # Very unbalanced cell types
           ref_matrix <- Matrix::Matrix(rpois(100 * 50, 10), 
                                        nrow = 100, ncol = 50, sparse = TRUE)
           query_matrix <- Matrix::Matrix(rpois(100 * 30, 10), 
                                          nrow = 100, ncol = 30, sparse = TRUE)
           
           # 80% Type1, 20% Type2
           ref_celltypes <- c(rep("Type1", 40), rep("Type2", 10))
           names(ref_celltypes) <- paste0("RefCell_", seq_len(50))
           
           query_celltypes <- c(rep("Type1", 24), rep("Type2", 6))
           names(query_celltypes) <- paste0("QueryCell_", seq_len(30))
           
           rownames(ref_matrix) <- paste0("Gene_", seq_len(100))
           colnames(ref_matrix) <- names(ref_celltypes)
           rownames(query_matrix) <- paste0("Gene_", seq_len(100))
           colnames(query_matrix) <- names(query_celltypes)
           
           list(
             ref_matrix = ref_matrix,
             ref_celltypes = ref_celltypes,
             query_matrix = query_matrix,
             query_celltypes = query_celltypes
           )
         },
         "sparse" = {
           # Very sparse data (lots of zeros)
           data <- create_small_test_data()
           # Make it more sparse by setting 70% to zero
           data$ref_matrix[sample(length(data$ref_matrix), 
                                  floor(0.7 * length(data$ref_matrix)))] <- 0
           data$query_matrix[sample(length(data$query_matrix), 
                                    floor(0.7 * length(data$query_matrix)))] <- 0
           data
         }
  )
}

#' Skip tests if running on CRAN or in limited environments
skip_if_no_parallel <- function() {
  if (!requireNamespace("parallel", quietly = TRUE)) {
    skip("parallel package not available")
  }
}

skip_if_slow <- function() {
  if (identical(Sys.getenv("NOT_CRAN"), "false")) {
    skip("Skipping slow test on CRAN")
  }
}