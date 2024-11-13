#' Run Constrained Supervised Factorization NMF
#'
#' @description
#' Performs constrained supervised factorization NMF analysis on the provided data
#'
#' @param object CSFNMF object
#' @param k Rank (number of factors). If NULL, estimated automatically
#' @param max_p_value Maximum p-value for rank selection
#' @param max.iter Maximum number of iterations
#' @param W0_H0 Initial W and H matrices
#' @param init_method Initialization method ("uniform", "regulated", "NNDSVD", "skmeanGenes", "skmeanCells")
#' @param theta Convergence threshold
#' @param const.alpha Alpha constraint parameter
#' @param const.beta Beta constraint parameter (NULL for automatic)
#' @param const.gamma Gamma constraint parameter
#' @param const.delta Delta constraint parameter
#' @param verbose Show progress messages
#' @param seed Random seed
#' @param num_cores Number of cores for parallel processing
#' @param chunk_size Size of chunks for processing (NULL for automatic)
#'
#' @return Updated CSFNMF object
#' @importFrom parallel mclapply detectCores
#' @importFrom progress progress_bar
#' @importFrom nnls nnls
#' @export
RunCSFNMF <- function(object,
                      k = NULL,
                      max_p_value = 0.01,
                      max.iter = 100,
                      W0_H0 = NULL,
                      init_method = "NNDSVD",
                      theta = 1.0e-8,
                      const.alpha = 1,
                      const.beta = NULL,
                      const.gamma = 1,
                      const.delta = 1,
                      verbose = TRUE,
                      seed = 1,
                      num_cores = 1,
                      chunk_size = NULL) {
  
  # Set up reporting
  report <- create_reporter(verbose)
  # Initialize training object
  report("Creating training object")
  train_object <- divide_reference_data(object, seed)
  
  # Select or set rank
  report("Determining rank")
  if (is.null(k)) {
    train_object@rank <- SelectRank(
      train_matrix = train_object@count.matrices@data,
      max_p_value = max_p_value,
      main_matrix = train_object@count.matrices@ref
    )
    k <- train_object@rank@k
  } else {
    train_object@rank@k <- k
    # Check k for skmeanGenes method
    if (init_method == "skmeanGenes") {
      max_k <- floor(nrow(train_object@count.matrices@ref) / 2) - 1
      if (k > max_k) {
        warning(sprintf("k reduced from %d to %d", k, max_k))
        train_object@rank@k <- max_k
      }
    }
  }
  
  
  # Set iterations
  train_object@max.iter <- max.iter
  
  # Initialize W and H if not provided
  if (is.null(W0_H0)) {
    report("Initializing W and H matrices")
    W0_H0 <- initialize_wh(train_object, init_method)
  }
  
  # Set matrices
  train_object@W <- W0_H0$W
  train_object@H <- W0_H0$H
  train_object@init_method <- init_method
  
  # Set beta constraint
  if (is.null(const.beta)) {
    const.beta <- standard_beta(train_object)
  } else if (is.numeric(const.beta)) {
    const.beta <- standard_beta(train_object, const.beta)
  }
  
  # Calculate helper matrices
  report("Calculating helper matrices")
  train_object <- calculate_help_matrices(train_object)
  
  # Calculate alpha
  report("Calculating alpha")
  real_alpha <- const.alpha
  real_N <- train_object@constants@N
  if (!is.numeric(const.alpha)) {
    train_object@constants@N <- const.alpha * train_object@constants@N
    const.alpha <- 1
  }
  
  # Set parameters
  train_object@hyper_para <- list(
    alpha = const.alpha,
    beta = const.beta,
    gamma = const.gamma,
    delta = const.delta
  )
  
  # Calculate constants for H
  report("Calculating H constants")
  train_object@constants@Hconst <- calculate_const_for_h(train_object)

  # Update W and H
  report("Updating W and H matrices")
  train_object <- update_wh(train_object, theta, verbose)
  
  # Restore original parameters
  train_object@hyper_para$alpha <- real_alpha
  train_object@constants@N <- real_N
  
  # Determine chunk size if not provided
  if (is.null(chunk_size)) {
    chunk_size <- determine_chunk_size(ncol(train_object@count.matrices@data))
  }
  
  # Calculate projection with optimized function
  report("Calculating H projection")
  h_project <- project_data(
    W = train_object@W,
    X = train_object@count.matrices@data,
    seed = seed,
    num_cores = num_cores,
    chunk_size = chunk_size,
    verbose = verbose
  )
  
  # Calculate accuracy
  report("Calculating accuracy")
  train_object@accuracy <- calculate_accuracy(train_object, h_project)
  if (verbose) {
    message(sprintf("Final accuracy: %.4f", train_object@accuracy))
  }
  
  # Update main object
  object@train_object <- train_object
  object@W <- ensure_dgCMatrix(train_object@W)
  object@H <- ensure_dgCMatrix(cbind(
    train_object@H,
    h_project
  )[, colnames(object@count.matrices@ref)])
  object@train_object <- train_object
  
  # Add processing info
  attr(object, "processing_info") <- list(
    num_cores_used = num_cores,
    chunk_size = chunk_size,
    iterations = length(train_object@loss),
    final_loss = tail(train_object@loss, 1),
    projection_info = attr(h_project, "processing_info")
  )
  
  object
}

#' Determine optimal chunk size
#'
#' @param n_cells Number of cells
#' @param available_memory Available memory in MB
#' @return Optimal chunk size
#' @keywords internal
determine_chunk_size <- function(n_cells, available_memory = 1000) {
  # Estimate memory per cell (adjust based on your data)
  mem_per_cell <- 0.1  # MB
  
  # Calculate optimal chunk size
  chunk_size <- min(
    floor(available_memory / mem_per_cell),
    n_cells,
    1000  # Maximum chunk size
  )
  
  # Ensure minimum size
  max(chunk_size, 100)
}

#' Ensure Matrix Format is dgCMatrix
#'
#' @param mat Input matrix
#' @return Matrix in dgCMatrix format
#' @keywords internal
ensure_dgCMatrix <- function(mat) {
  if (!inherits(mat, "dgCMatrix")) {
    return(as(as(mat, "CsparseMatrix"), "dgCMatrix"))
  }
  return(mat)
}