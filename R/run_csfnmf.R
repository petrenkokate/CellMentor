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
#' @param init_method Initialization method
#' @param theta Convergence threshold
#' @param const.alpha Alpha constraint parameter
#' @param const.beta Beta constraint parameter
#' @param const.gamma Gamma constraint parameter
#' @param const.delta Delta constraint parameter
#' @param verbose Show progress messages
#' @param seed Random seed
#' @param num_cores Number of cores for parallel processing
#'
#' @return Updated CSFNMF object
#' @export
RunCSFNMF <- function(train_object,
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
                      whole_object = FALSE) {
  
  report <- create_reporter(verbose)
  
  
  # Initialize parameters
  parameters <- list(
    rank = k,
    max_iter = max.iter,
    init_method = init_method,
    alpha = const.alpha,
    beta = const.beta,
    gamma = const.gamma,
    delta = const.delta
  )
  train_object@parameters <- parameters
  # Initialize results list
  train_object@results <- list(
    loss = numeric(0),
    nmi = NA,
    predictions = NULL
  )

  # Initialize W and H if not provided
  if (is.null(W0_H0)) {
    report("Initializing W and H matrices")
    W0_H0 <- initialize_wh(train_object, init_method)
  }
  
  # Set matrices
  train_object@W <- W0_H0$W
  train_object@H <- W0_H0$H
  
  # Set beta constraint
  if (is.null(const.beta)) {
    parameters$beta <- standard_beta(train_object)
  } else if (is.numeric(const.beta)) {
    parameters$beta <- standard_beta(train_object, const.beta)
  }
  
  report("Calculating helper matrices")
  train_object@constants <- calculate_help_matrices(train_object)
  
  # Calculate alpha
  report("Calculating alpha")
  real_alpha <- const.alpha
  real_N <- train_object@constants@N
  if (!is.numeric(const.alpha)) {
    train_object@constants@N <- const.alpha * train_object@constants@N
    parameters$alpha <- 1
  }
  
  # Set parameters
  train_object@parameters <- parameters
  
  # Calculate constants for H
  report("Calculating H constants")
  train_object@constants@Hconst <- calculate_const_for_h(train_object)
  
  # Update W and H
  report("Updating W and H matrices")
  train_object <- update_wh(train_object, theta, verbose)
  
  # Restore original parameters
  parameters$alpha <- real_alpha
  train_object@constants@N <- real_N
  
  if (!whole_object) {
    # Project data and calculate performance
    report("Calculating projections and performance metrics")
    h_project <- project_data(
      W = train_object@W,
      X = train_object@matrices@data,
      seed = seed,
      num_cores = num_cores,
      verbose = verbose
    )
    
    # Calculate NMI
    report("Calculating NMI")
    perf_result <- calculate_performance(train_object, h_project, return_pred = TRUE)
    train_object@results$nmi <- perf_result$nmi
    train_object@results$predictions <- perf_result$SingleRpred
    
    if (verbose) {
      message(sprintf("Final NMI: %.4f", perf_result$nmi))
    }
    # Add processing info
    attr(train_object, "processing_info") <- list(
      num_cores_used = num_cores,
      iterations = length(train_object@results$loss),
      final_loss = tail(train_object@results$loss, 1),
      final_nmi = perf_result$nmi,
      projection_info = attr(h_project, "processing_info")
    )
    
  }
  
  else {
    # Add processing info
    attr(train_object, "processing_info") <- list(
      num_cores_used = num_cores,
      iterations = length(train_object@results$loss),
      final_loss = tail(train_object@results$loss, 1)
    )
  }
  
  train_object
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