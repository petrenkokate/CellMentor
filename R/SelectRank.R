#' @importFrom methods setClass new
#' @importFrom stats ks.test median rnorm
#' @importFrom RMTstat qmp dmp
#' @importFrom Matrix Matrix t crossprod tcrossprod Diagonal
#' @importFrom parallel mclapply
#' @importFrom ggplot2 ggplot aes geom_histogram ggplot_build
#' @importFrom data.table as.data.table
#' @importFrom magrittr %>%
NULL


#' Select Optimal Rank
#'
#' @param train_matrix Training matrix
#' @param max_p_value Maximum p-value threshold (default: 0.01)
#' @param jump_step Step size for beta search (default: 0.05)
#' @param numCores Number of cores for parallel processing (default: 1)
#' @param range Range for beta search (default: c(0,1))
#' @param main_matrix Optional main matrix
#' @param verbose Show progress messages
#' @return ExtendedRank object
#' @export
SelectRank <- function(train_matrix,
                       max_p_value = 0.01,
                       jump_step = 0.05,
                       numCores = 1,
                       range = c(0, 1),
                       main_matrix = NULL,
                       verbose = TRUE) {
  
  report <- create_reporter(verbose)
  main_matrix <- main_matrix %||% train_matrix
  
  # Calculate optimal beta values
  report("Calculating optimal beta values")
  beta_list <- seq(range[1], range[2], by = jump_step)
  
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
  
  # Initialize parameters
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
  
  methods::new(
    Class = "ExtendedRank",
    k = count_k,
    alpha = alpha,
    beta = beta,
    scale_vectors = scales
  )
}

#' Calculate Significant Components
#' @keywords internal
calculate_significant_components <- function(vectors, k, n.genes, max_p_value) {
  sum(vapply(seq_len(k), function(i) {
    vector <- vectors[,i]
    p_value <- stats::ks.test(
      vector,
      rnorm(length(vector), mean = 0, sd = sqrt(1/n.genes))
    )[["p.value"]]
    p_value <= max_p_value
  }, logical(1)))
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


#' Sinkhorn-Knopp Matrix Scaling
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
  
  # Prepare matrix
  var_A <- if (is.null(alpha) || is.null(beta)) {
    A
  } else {
    alpha * ((1-beta) * A + beta * A^2)
  }
  
  # Initialize
  n <- ncol(A)
  m <- nrow(A)
  x <- matrix(1, nrow = m, ncol = 1)
  y <- matrix(1, nrow = n, ncol = 1)
  
  # Iterate
  for (tau in seq_len(niter)) {
    # Update vectors with stabilization
    y_new <- stabilize_division(sum_col, Matrix::t(var_A) %*% x)
    x_new <- stabilize_division(sum_row, var_A %*% y_new)
    
    # Check convergence
    converged <- check_convergence(
      x_new, y_new, var_A, sum_row, sum_col, tol
    )
    
    if (converged) {
      return(list(
        x = pmax(as.numeric(x_new), 1e-10),
        y = pmax(as.numeric(y_new), 1e-10)
      ))
    }
    
    x <- x_new
    y <- y_new
  }
  
  list(x = as.numeric(x), y = as.numeric(y))
}

#' Select Alpha for Constant Beta
#'
#' @param Y Input matrix
#' @param beta Beta value
#' @return List with alpha and scale matrix
#' @keywords internal
SelectAlphaForConstBeta <- function(Y, beta) {
  Y <- ensure_matrix(Y)
  num_cells <- ncol(Y)
  num_genes <- nrow(Y)
  
  u_v <- NonSquareSinkhornKnopp(Y, alpha = 1, beta = beta)
  
  x_vec <- pmax(u_v[["x"]], 1e-10)
  y_vec <- pmax(u_v[["y"]], 1e-10)
  
  Y_hat <- with_memory_check({
    diag(sqrt(x_vec), nrow = num_genes, ncol = num_genes) %*% Y %*% 
      diag(sqrt(y_vec), nrow = num_cells, ncol = num_cells)
  })

  dubble_Y_hat <- with_memory_check({
    (1/num_cells) * (Y_hat %*% Matrix::t(Y_hat) + diag(1e-10, num_genes))
  })

  ev_list <- eigen(dubble_Y_hat)$values
  gamma <- num_genes/num_cells
  
  lambda_median <- max(median(ev_list), 1e-10)
  meu_half <- max(RMTstat::qmp(p = 1/2, svr = 1/gamma, var = 1), 1e-10)
  
  list(
    alpha = lambda_median/meu_half,
    scale_matrix = dubble_Y_hat
  )
}


#' Find Optimal Beta
#'
#' @param data Input matrix
#' @param beta Beta value
#' @return Difference statistic
#' @keywords internal
OptimalBeta <- function(data, beta) {
  data <- ensure_matrix(data)
  num_cells <- ncol(data)
  num_genes <- nrow(data)
  
  result <- SelectAlphaForConstBeta(data, beta)
  dubble_data_hat <- (1/result$alpha) * result$scale_matrix
  
  ev_hat <- as.data.frame(eigen(dubble_data_hat)$values)
  colnames(ev_hat) <- "value"
  
  p <-  ggplot2::ggplot(ev_hat, aes(x = value)) +
    ggplot2::geom_histogram(aes(y = ..density..), binwidth = 0.05)
  
  x <- data.table::as.data.table((ggplot2::ggplot_build(p)$data[1]))
  mp <- RMTstat::dmp(
    x = x$x,
    svr = num_cells/num_genes,
    var = 1
  )
  
  stats::ks.test(x$y, mp)[["statistic"]][["D"]]
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
  data_matrix <- object@count.matrices@ref
  
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