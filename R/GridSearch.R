#' Select Optimal Hyperparameters for CSFNMF
#'
#' @description
#' Tests different combinations of hyperparameters using `RunCSFNMF` to find the optimal configuration.
#' The function performs a grid search over specified parameter ranges, evaluating the model's 
#' performance for each combination. Parameters `alpha` and `beta` are kept equal during optimization.
#' The rank (`k`) can be provided, taken from an existing object, or determined automatically using
#' `SelectRank`.
#'
#' @param object CSFNMF object containing reference and query data matrices, with required matrices 
#'               for `data` and `ref` under `object@matrices`.
#' @param k Optional rank value (number of factors). If `NULL`:
#'          - Uses existing rank from `object` if available
#'          - Otherwise determines automatically using `SelectRank`
#' @param init_methods Vector of initialization methods to test. Options:
#'          - "uniform": Random uniform initialization
#'          - "regulated": Cell-type guided initialization
#'          - "NNDSVD": Non-negative Double SVD
#'          - "skmeanGenes": Gene clustering-based
#'          - "skmeanCells": Cell clustering-based
#'          Default: c("NNDSVD")
#' @param alpha_range Vector of `alpha` values to test.
#'          Controls within-class scatter (cell similarity within the same type).
#'          Default: c(1, 5, 10)
#' @param beta_range Vector of `beta` values to test.
#'          Controls between-class scatter (cell separation between different types).
#'          Default: c(1, 5, 10)
#' @param gamma_range Vector of sparsity parameter values to test.
#'          Controls sparsity of the factorization.
#'          Default: c(0, 0.001, 0.1, 1)
#' @param delta_range Vector of orthogonality parameter values to test.
#'          Controls orthogonality between factors.
#'          Default: c(0, 0.1, 1, 2.5)
#' @param n_iter Number of repetitions per configuration for averaging results (default: 3).
#' @param verbose Logical; whether to show progress messages during optimization.
#'          Default: `TRUE`
#' @param num_cores Number of cores to use for parallel processing.
#'          If > 1, parameter combinations are tested in parallel.
#'          Default: 1
#' @param subset_size Optional integer specifying the number of cells to use in a subset for faster
#'          parameter optimization. When provided, a random subset of cells is selected from the data
#'          and used to test different parameter configurations. This can reduce computation time,
#'          but the best parameters found may vary slightly from those found with the full dataset.
#'          After selecting the optimal parameters, the function runs the final model on the full
#'          dataset. Default: `NULL` (uses full dataset).
#'
#' @return List containing:
#'         - `best_params`: List with the overall best parameter configuration:
#'           * `k`: Selected rank
#'           * `init_method`: Best initialization method
#'           * `alpha`: Best alpha parameter
#'           * `beta`: Best beta parameter
#'           * `gamma`: Best gamma value
#'           * `delta`: Best delta value
#'           * `accuracy`: Best achieved accuracy
#'           * `loss`: Corresponding loss value
#'         - `results`: Data frame of all combinations tested, including:
#'           * `init_method`: Initialization method used
#'           * `alpha`: Alpha parameter value
#'           * `beta`: Beta parameter value
#'           * `gamma`: Gamma parameter value
#'           * `delta`: Delta parameter value
#'           * `accuracy`: Achieved accuracy
#'           * `loss`: Final loss value
#'           * `convergence_iter`: Number of iterations for convergence
#'         - `best_model`: CSFNMF model object trained with the best parameters.
#'
#' @examples
#' # Basic usage with default parameter ranges
#' optimal_params <- select_optimal_parameters(object)
#'
#' # Custom parameter ranges with parallel processing
#' optimal_params <- select_optimal_parameters(
#'   object,
#'   k = 50,
#'   init_methods = c("NNDSVD"),
#'   alpha_range = c(1, 2, 3),
#'   beta_range = c(1, 2, 5),
#'   gamma_range = c(0, 0.1, 0.2),
#'   delta_range = c(0, 1),
#'   num_cores = 4
#' )
#'
#' # Access best overall parameters
#' print(optimal_params$best_params)
#'
#' # Compare results across methods
#' lapply(optimal_params$results, function(x) x$best_params)
#'
#' @seealso
#' \code{\link{RunCSFNMF}} for the main model fitting function
#' \code{\link{SelectRank}} for automatic rank selection
#' \code{\link{plot_parameter_search}} for results visualization
#' 
#' @export
select_optimal_parameters <- function(object,
                                      k = NULL,
                                      init_methods = c("uniform", "regulated", "NNDSVD", 
                                                       "skmeanGenes", "skmeanCells"),
                                      alpha_range = c(1, 5, 10),
                                      beta_range = c(1, 5, 10),
                                      gamma_range = c(0, 0.001, 0.1, 1),
                                      delta_range = c(0, 0.1, 1, 2.5),
                                      n_iter = 1,
                                      verbose = TRUE,
                                      num_cores = 1,
                                      subset_size = 0.3) {
  
  # Create reporter for progress messages
  report <- create_reporter(verbose)
  
  # Handle data subsetting if requested
  working_object <- if(subset_size < 1) {
    report(sprintf("Creating %.0f%% data subset for parameter search", subset_size * 100))
    create_data_subset(object, subset_size)
  } else {
    report("Using full dataset for parameter search")
    object
  }
  
  # Determine k value
  if (is.null(k)) {
    if (!is.null(object@train_object) && !is.null(object@train_object@parameters$rank) &&
        object@train_object@parameters$rank > 0) {
      k <- object@train_object@parameters$rank
      report(sprintf("Using existing rank: %d", k))
    } else {
      report("Determining optimal rank")
      rank_result <- SelectRank(
        train_matrix = working_object@matrices@data,
        max_p_value = 0.01,
        main_matrix = working_object@matrices@ref,
        verbose = verbose,
        numCores = num_cores
      )
      k <- rank_result$rank
      report(sprintf("Optimal rank determined: %d", k))
    }
  } else {
    report(sprintf("Using provided rank: %d", k))
  }
  
  # Create parameter grid
  param_grid <- expand.grid(
    init_method = init_methods,
    alpha = alpha_range,
    beta = beta_range,
    gamma = gamma_range,
    delta = delta_range,
    stringsAsFactors = FALSE
  )
  
  # Initialize results storage
  results <- data.frame(
    init_method = character(),
    alpha = numeric(),
    beta = numeric(),
    gamma = numeric(),
    delta = numeric(),
    accuracy = numeric(),
    loss = numeric(),
    convergence_iter = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Function to evaluate one parameter combination
  evaluate_params <- function(params) {
    tryCatch({
      model <- RunCSFNMF(
        object = working_object,
        k = k,
        init_method = params$init_method,
        const.alpha = params$alpha,
        const.beta = params$beta,
        const.gamma = params$gamma,
        const.delta = params$delta,
        max.iter = 100,
        verbose = FALSE,
        num_cores = num_cores
      )
      
      c(
        accuracy = model@train_object@results$accuracy,
        loss = tail(model@train_object@results$loss, 1),
        convergence_iter = length(model@train_object@results$loss)
      )
    }, error = function(e) {
      warning("Error with parameters: ", paste(params, collapse = ", "))
      c(accuracy = NA, loss = NA, convergence_iter = NA)
    })
  }
  
  # Evaluate all parameter combinations
  report("Starting parameter grid search")
  total_configs <- nrow(param_grid)
  
  for (i in seq_len(total_configs)) {
    if (verbose) {
      report(sprintf("Testing configuration %d/%d", i, total_configs))
      print(param_grid[i, ])
    }
    
    iter_results <- replicate(n_iter, {
      res <- evaluate_params(param_grid[i, ])
      res[1:3]
    }, simplify = TRUE)
    
    if (n_iter == 1) {
      mean_results <- iter_results
    } else {
      mean_results <- rowMeans(iter_results, na.rm = TRUE)
    }
    
    results <- rbind(results, data.frame(
      param_grid[i, ],
      accuracy = mean_results["accuracy",],
      loss = mean_results["loss",],
      convergence_iter = mean_results["convergence_iter",]
    ))
  }
  
  # Find best parameters
  best_idx <- which.max(results$accuracy)
  best_params <- param_grid[best_idx, ]
  
  # Run final model with best parameters on full dataset
  report("Training final model with best parameters on full dataset")
  best_model <- RunCSFNMF(
      object = object,
      k = k,
      init_method = best_params$init_method,
      const.alpha = best_params$alpha,
      const.beta = best_params$beta,
      const.gamma = best_params$gamma,
      const.delta = best_params$delta,
      max.iter = 100,
      verbose = verbose,
      num_cores = num_cores
    )
  
  list(
    best_params = c(k = k, best_params),
    results = results,
    best_model = best_model
  )
}

#' Helper function to create data subset
#' 
#' @export
create_data_subset <- function(object, subset_size) {
  # Calculate number of cells to sample
  n_ref_cells <- ncol(object@matrices@ref)
  n_data_cells <- ncol(object@matrices@data)
  n_ref_subset <- floor(n_ref_cells * subset_size)
  n_data_subset <- floor(n_data_cells * subset_size)
  
  # Sample cells
  ref_cells <- sample(seq_len(n_ref_cells), n_ref_subset)
  data_cells <- sample(seq_len(n_data_cells), n_data_subset)
  
  # Create new object with subset
  new_object <- methods::new("csfnmf")
  new_object@matrices <- methods::new(
    "RefDataList",
    ref = object@matrices@ref[, ref_cells],
    data = object@matrices@data[, data_cells]
  )
  new_object@annotation <- object@annotation[ref_cells, , drop = FALSE]
  
  new_object
}

#' Plot Parameter Search Results
#'
#' @param results Results from select_optimal_parameters
#' @return ggplot object with visualization of parameter search results
#' @export
plot_parameter_search <- function(results) {
  # Create accuracy heatmap for each initialization method
  ggplot2::ggplot(results$results, 
                  aes(x = factor(alpha), y = factor(beta), fill = accuracy)) +
    geom_tile() +
    facet_wrap(~init_method) +
    scale_fill_viridis_c() +
    labs(title = "Parameter Search Results",
         x = "Alpha",
         y = "Beta",
         fill = "Accuracy") +
    theme_minimal()
}