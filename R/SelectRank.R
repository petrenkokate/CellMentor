#' @importFrom methods setClass new
#' @importFrom stats ks.test median rnorm
#' @importFrom RMTstat qmp dmp
#' @importFrom Matrix Matrix t crossprod tcrossprod Diagonal
#' @importFrom parallel mclapply
#' @importFrom ggplot2 ggplot aes geom_histogram ggplot_build
#' @importFrom data.table as.data.table
#' @importFrom magrittr %>%
NULL

#' Optimized Rank Selection
#'
#' @param train_matrix Training matrix
#' @param max_p_value Maximum p-value threshold (default: 0.01)
#' @param jump_step Step size for beta search (default: 0.05)
#' @param numCores Number of cores for parallel processing (default: 1)
#' @param range Range for beta search (default: c(0,1))
#' @param main_matrix Optional main matrix
#' @param verbose Show progress
#' @return List containing rank value and scaling parameters
#' @importFrom Matrix Matrix t crossprod tcrossprod Diagonal
#' @export
SelectRank <- function(train_matrix,
                                 max_p_value = 0.01,
                                 jump_step = 0.05,
                                 numCores = 1,
                                 range = c(0, 1),
                                 verbose = TRUE) {
  
  report <- if(verbose) message else function(...) NULL
  main_matrix <- train_matrix
  
  # Pre-compute matrix dimensions and constants
  n.genes <- nrow(main_matrix)
  n.cells <- ncol(main_matrix)
  upper_edge <- (1 + sqrt(n.genes/n.cells))^2
  
  # Optimize beta search
  beta_list <- seq(range[1], range[2], by = jump_step)
  
  # Pre-allocate results storage
  best_params <- list(k = 0, alpha = 0, beta = 0, scales = NULL)
  
  diff_list <- if (numCores > 1) {
    parallel::mclapply(
      beta_list,
      function(x) OptimalBeta(train_matrix, x),
      mc.cores = min(length(beta_list), numCores)
    )
  } else {
    lapply(beta_list, function(x) OptimalBeta(train_matrix, x))
  }
  
  diff_values <- unlist(diff_list)
  beta_opt_list <- beta_list[which(diff_values == min(diff_values))]
  
  k <- beta <- alpha <- 0
  scales <- list()
  n.genes <- nrow(main_matrix)
  n.cells <- ncol(main_matrix)
  upper_edge <- (1 + sqrt(n.genes/n.cells))^2
  
  report(sprintf("Upper edge: %.2f", upper_edge))
  
  # Process optimal beta values
  for (beta_opt in beta_opt_list) {
    result <- with_memory_check({
      SelectAlphaForConstBeta(train_matrix, beta = beta_opt)
    })
    
    u_v <- with_memory_check({
      NonSquareSinkhornKnopp(main_matrix, alpha = result$alpha, beta = beta_opt)
    })
    
    # Stabilize values and calculate matrices
    x_vec <- pmax(u_v[["x"]], 1e-10)
    y_vec <- pmax(u_v[["y"]], 1e-10)
    
    X_hat <- with_memory_check({
      diag(sqrt(x_vec), nrow = n.genes, ncol = n.genes) %*%
        main_matrix %*%
        diag(sqrt(y_vec), nrow = n.cells, ncol = n.cells)
    })
    
    dubble_X_hat <- with_memory_check({
      (1/n.cells) * (X_hat %*% Matrix::t(X_hat) + diag(1e-10, n.genes))
    })
    
    ev_hat <- eigen(dubble_X_hat)
    k_opt <- sum(ev_hat$values > upper_edge)
    
    report(sprintf("Initial k_opt based on eigenvalues: %d", k_opt))
    
    if (k_opt > k) {
      k <- k_opt
      alpha <- result$alpha
      beta <- beta_opt
      scales <- u_v
    }
  }
  
  
  # Calculate final k
  count_k <- calculate_significant_components(
    ev_hat$vectors,
    k,
    n.genes,
    max_p_value
  )
  
  report(sprintf("Final k after p-value filtering: %d", count_k))
  
  count_k
}

#' Optimized Component Significance Calculation
#' @keywords internal
calculate_significant_components <- function(vectors, k, n.genes, max_p_value) {
  # Pre-generate normal distribution parameters
  norm_mean <- 0
  norm_sd <- sqrt(1/n.genes)
  
  # Vectorized computation
  p_values <- vapply(seq_len(k), function(i) {
    stats::ks.test(
      vectors[,i],
      "pnorm",
      mean = norm_mean,
      sd = norm_sd
    )$p.value
  }, numeric(1))
  
  sum(p_values <= max_p_value)
}

#' Ensure Matrix Format
#' 
#' @param mat Input matrix
#' @return Properly formatted sparse matrix
#' @importFrom Matrix sparse.model.matrix
#' @keywords internal
ensure_matrix <- function(mat) {
  tryCatch({
    if (is(mat, "Matrix")) {
      as(mat, "CsparseMatrix")
    } else if (is.matrix(mat)) {
      as(mat, "CsparseMatrix")
    } else {
      as(as.matrix(mat), "CsparseMatrix")
    }
  }, error = function(e) {
    stop("Could not convert input to sparse matrix format: ", e$message)
  })
}

#' Print Matrix Info
#' 
#' @param mat Matrix to examine
#' @param name Name for output
#' @keywords internal
debug_matrix <- function(mat, name = "Matrix") {
  message(sprintf(
    "%s: [%d x %d] %s, %.2f%% sparse",
    name,
    nrow(mat),
    ncol(mat),
    class(mat)[1],
    100 * (1 - length(mat@x)/(nrow(mat) * ncol(mat)))
  ))
}

#' Corrected Sinkhorn-Knopp Matrix Scaling
#'
#' @param A Input matrix
#' @param sum_row Row sums
#' @param sum_col Column sums
#' @param niter Maximum iterations
#' @param tol Convergence tolerance
#' @param alpha Alpha parameter
#' @param beta Beta parameter
#' @return List of scaling vectors
#' @importFrom Matrix t crossprod tcrossprod Diagonal
#' @keywords internal
NonSquareSinkhornKnopp <- function(A,
                                             sum_row = matrix(ncol(A), nrow = nrow(A), ncol = 1),
                                             sum_col = matrix(nrow(A), nrow = ncol(A), ncol = 1),
                                             niter = 100,
                                             tol = 1e-10,
                                             alpha = NULL,
                                             beta = NULL) {
  
  # Matrix preparation with corrected quadratic term
  var_A <- if (is.null(alpha) || is.null(beta)) {
    A
  } else {
    # For element-wise operations, use direct multiplication
    A_squared <- A * A  # Element-wise multiplication
    alpha * ((1-beta) * A + beta * A_squared)
  }
  
  # Pre-compute dimensions
  n <- ncol(A)
  m <- nrow(A)
  
  # Initialize scaling vectors
  x <- rep(1, m)
  y <- rep(1, n)
  
  # Pre-compute transpose
  var_A_t <- Matrix::t(var_A)
  
  # Convergence tracking
  conv_history <- numeric(3)
  eps <- 1e-10  # Stabilization factor
  
  # Optimized iteration
  for (tau in seq_len(niter)) {
    # Update vectors with stabilization
    y_new <- as.numeric(sum_col/(var_A_t %*% x + eps))
    x_new <- as.numeric(sum_row/(var_A %*% y_new + eps))
    
    # Check convergence efficiently
    x_diff <- abs(x_new - x)
    y_diff <- abs(y_new - y)
    max_diff <- max(c(x_diff, y_diff))
    
    # Update convergence history
    conv_history <- c(conv_history[-1], max_diff)
    
    # Check convergence with acceleration
    if (max_diff < tol || 
        (tau > 3 && all(diff(conv_history) > -tol/10))) {
      return(list(
        x = pmax(x_new, eps),
        y = pmax(y_new, eps)
      ))
    }
    
    x <- x_new
    y <- y_new
  }
  
  list(x = pmax(x, eps), y = pmax(y, eps))
}


#' Optimized Alpha Selection for Constant Beta
#'
#' @param Y Input matrix
#' @param beta Beta value
#' @return List with alpha and scale matrix
#' @importFrom Matrix t crossprod tcrossprod Diagonal
#' @keywords internal
SelectAlphaForConstBeta <- function(Y, beta) {
  # Initial matrix preparation
  Y <- ensure_matrix(Y)
  num_cells <- ncol(Y)
  num_genes <- nrow(Y)
  
  # Get scaling vectors using optimized Sinkhorn-Knopp
  u_v <- NonSquareSinkhornKnopp(Y, alpha = 1, beta = beta)
  
  # Stabilize scaling vectors
  x_vec <- pmax(u_v[["x"]], 1e-10)
  y_vec <- pmax(u_v[["y"]], 1e-10)
  
  # Create diagonal scaling matrices
  Dx <- Matrix::Diagonal(n = num_genes, x = sqrt(x_vec))
  Dy <- Matrix::Diagonal(n = num_cells, x = sqrt(y_vec))
  
  # Compute scaled matrix Y_hat
  Y_hat <- Dx %*% Y %*% Dy
  
  # Compute covariance matrix with stabilization
  dubble_Y_hat <- (1/num_cells) * (
    Y_hat %*% Matrix::t(Y_hat) + 
      Matrix::Diagonal(n = num_genes, x = 1e-10)
  )
  
  # Ensure symmetric property for eigenvalue computation
  dubble_Y_hat <- (dubble_Y_hat + Matrix::t(dubble_Y_hat))/2
  
  # Compute eigenvalues
  ev_list <- eigen(dubble_Y_hat, symmetric = TRUE)$values
  
  # Calculate scaling parameters
  gamma <- num_genes/num_cells
  lambda_median <- max(median(ev_list), 1e-10)
  meu_half <- max(RMTstat::qmp(p = 1/2, svr = 1/gamma, var = 1), 1e-10)
  
  # Return results
  list(
    alpha = lambda_median/meu_half,
    scale_matrix = dubble_Y_hat
  )
}


#' Optimized Optimal Beta Calculation
#'
#' @param data Input matrix
#' @param beta Beta value
#' @return Difference statistic
#' @importFrom stats hist
#' @keywords internal
OptimalBeta <- function(data, beta) {
  # Initial setup
  data <- ensure_matrix(data)
  num_cells <- ncol(data)
  num_genes <- nrow(data)
  
  # Get eigenvalues
  result <- SelectAlphaForConstBeta(data, beta)
  dubble_data_hat <- (1/result$alpha) * result$scale_matrix
  eigenvalues <- eigen(dubble_data_hat)$values
  
  # Create fixed-width bins exactly as ggplot2 does
  range_ev <- range(eigenvalues)
  binwidth <- 0.05
  bins <- seq(floor(range_ev[1]/binwidth) * binwidth,
              ceiling(range_ev[2]/binwidth) * binwidth,
              by = binwidth)
  
  # Calculate histogram
  h <- hist(eigenvalues, breaks = bins, plot = FALSE)
  
  # Get density values at bin centers (equivalent to ggplot2's calculation)
  x_values <- h$mids
  y_values <- h$density
  
  # Calculate MP density at the same points
  mp_density <- RMTstat::dmp(
    x = x_values,
    svr = num_cells/num_genes,
    var = 1
  )
  
  # Perform KS test
  ks_result <- stats::ks.test(y_values, mp_density)
  
  return(ks_result$statistic[["D"]])
}

#' Helper Functions
#' @keywords internal
stabilize_division <- function(num, denom, eps = 1e-10) {
  num/(denom + eps)
}

#' @keywords internal
check_convergence <- function(x_new, y_new, var_A, sum_row, sum_col, tol) {
  m <- nrow(var_A)
  n <- ncol(var_A)
  
  max_i <- max(abs(diag(as.numeric(x_new), nrow = m, ncol = m) %*%
                     var_A %*% y_new - sum_row), na.rm = TRUE)
  max_j <- max(abs(Matrix::t(x_new) %*% var_A %*%
                     diag(as.numeric(y_new), nrow = n, ncol = n) - Matrix::t(sum_col)),
               na.rm = TRUE)
  
  max_i < tol && max_j < tol
}

#' @keywords internal
calculate_edges <- function(gamma) {
  list(
    upper = (1 + sqrt(gamma))^2,
    lower = (1 - sqrt(gamma))^2
  )
}

#' Matrix Connectivity Analysis
#' @export
ConnectivityGroups <- function(data_matrix) {
  genes_Q <- rownames(data_matrix)
  cells_Q <- colnames(data_matrix)
  connectivity <- data.frame(
    group_num = rep(0, length(cells_Q)),
    row.names = cells_Q
  )
  
  process_groups(data_matrix, genes_Q, cells_Q, connectivity)
}

#' @keywords internal
process_groups <- function(data_matrix, genes_Q, cells_Q, connectivity) {
  group <- 1
  
  while (length(cells_Q) > 0) {
    # Process current group
    result <- process_single_group(
      data_matrix, genes_Q, cells_Q, connectivity, group
    )
    
    connectivity <- result$connectivity
    cells_Q <- result$remaining_cells
    group <- group + 1
  }
  
  connectivity
}

#' Select Cells for Finding K
#' @export
SelectCellsForFindK <- function(object, gamma = 0.5, seed = 1) {
  set.seed(seed)
  # Change from count.matrices to matrices
  data_matrix <- object@matrices@ref
  
  # Calculate parameters
  num_cells <- ncol(data_matrix)
  new_num_cells <- nrow(data_matrix) * (1/gamma)
  cell_counts <- table(object@annotation$celltype)
  
  # Select cells by type
  selected_cells <- select_cells_by_type(
    data_matrix, object@annotation$celltype,
    cell_counts, new_num_cells, num_cells
  )
  
  # Handle zero expression genes
  selected_cells <- handle_zero_genes(data_matrix, selected_cells)
  
  data_matrix[, selected_cells]
}

#' @keywords internal
select_cells_by_type <- function(data_matrix, cell_types, counts, new_total, current_total) {
  unlist(lapply(names(counts), function(type) {
    cells <- names(which(cell_types == type))
    n_select <- ceiling((counts[type] * new_total)/current_total)
    sample(cells, n_select)
  }))
}

#' @keywords internal
handle_zero_genes <- function(data_matrix, selected_cells) {
  new_data <- data_matrix[, selected_cells]
  zero_genes <- rownames(new_data)[Matrix::rowSums(new_data) == 0]
  
  if (length(zero_genes) > 0) {
    additional_cells <- unlist(lapply(zero_genes, function(gene) {
      expressing_cells <- which(data_matrix[gene,] > 0)
      if (length(expressing_cells) > 0) {
        sample(colnames(data_matrix)[expressing_cells],
               min(5, length(expressing_cells)))
      }
    }))
    
    selected_cells <- unique(c(selected_cells, additional_cells))
  }
  
  selected_cells
}