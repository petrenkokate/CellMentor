#' Calculate Loss Function
#'
#' @param object CSFNMF object
#' @param divided Logical: return components separately
#' @return Loss value or list of components
#' @importFrom Matrix norm crossprod
#' @keywords internal
calculate_loss <- function(object, divided = FALSE) {
  # Extract matrices and parameters
  X <- object@matrices@ref
  W <- object@W
  H <- object@H
  const.list <- object@parameters
  num.clusters <- length(table(object@annotation$celltype))
  
  # Pre-compute scaling factors
  alpha_factor <- const.list[["alpha"]]/(2 * num.clusters)
  cluster_factor <- 1/(2 * (num.clusters - 1) * num.clusters)
  
  # Calculate reconstruction error using efficient matrix operations
  WH <- W %*% H
  error <- norm(X - WH, type = "F")^2
  
  # Calculate within-cluster variance efficiently
  H_diff <- H - H %*% object@constants@A
  within_var <- sum(H_diff * (H_diff %*% object@constants@N))
  
  # Calculate between-cluster variance
  H_between <- H %*% object@constants@BP
  between_var <- sum(H_between * (H_between %*% const.list[["beta"]]))
  
  # Calculate regularization terms
  reg_h <- sum(H)
  reg_w <- sum(crossprod(W) * object@constants@M)
  
  if (divided) {
    list(
      error = 0.5 * error,
      sw = alpha_factor * within_var,
      sb = cluster_factor * between_var,
      reg.h = const.list[["gamma"]] * reg_h,
      reg.w = (const.list[["delta"]]/2) * reg_w
    )
  } else {
    0.5 * error + 
      alpha_factor * within_var - 
      cluster_factor * between_var + 
      const.list[["gamma"]] * reg_h + 
      (const.list[["delta"]]/2) * reg_w
  }
}

#' Update W Matrix
#'
#' @param object CSFNMF object
#' @return Updated W matrix
#' @importFrom Matrix crossprod tcrossprod
#' @keywords internal
update_w <- function(object) {
  eps <- .Machine$double.eps
  
  # Optimize numerator computation
  # X %*% t(H) using efficient matrix multiplication
  up <- object@matrices@ref %*% t(object@H)
  
  # Optimize denominator computation
  # W %*% (H %*% t(H) + delta * M)
  HtH <- tcrossprod(object@H)
  reg_term <- object@parameters[["delta"]] * object@constants@M
  down <- object@W %*% (HtH + reg_term)
  
  # Element-wise operations for multiplicative update
  W_new <- object@W * (up/(down + eps))
  
  # Return as dense matrix
  as(W_new, "dgeMatrix")
}

#' Update H Matrix
#'
#' @param object CSFNMF object
#' @param num_cores Number of cores for parallel processing
#' @return Updated H matrix
#' @importFrom parallel mclapply
#' @keywords internal
update_h <- function(object, num_cores = 1) {
  eps <- .Machine$double.eps
  
  # Pre-compute constant matrices
  WtW <- crossprod(object@W, object@W)
  WtX <- crossprod(object@W, object@matrices@ref)
  gamma_matrix <- matrix(object@parameters[["gamma"]], 
                         nrow = nrow(object@H), 
                         ncol = ncol(object@H))
  
  # Compute core components efficiently
  H_pos <- object@H %*% object@constants@Hconst[["positive"]]
  H_neg <- object@H %*% object@constants@Hconst[["negative"]]
  
  if (num_cores > 1) {
    # Optimize parallel processing
    num_cores <- min(num_cores, parallel::detectCores(), ncol(object@H))
    chunk_size <- ceiling(ncol(object@H) / num_cores)
    chunks <- split(seq_len(ncol(object@H)), 
                    ceiling(seq_len(ncol(object@H))/chunk_size))
    
    # Process chunks in parallel with optimized operations
    H_new <- do.call(cbind, parallel::mclapply(chunks, function(idx) {
      # Efficient computation for each chunk
      a <- WtW %*% object@H[,idx] + H_pos[,idx] + gamma_matrix[,idx]
      b <- WtX[,idx]
      c <- H_neg[,idx]
      
      # Compute update using quadratic formula
      up <- b + sqrt(b * b + 4 * a * c)
      down <- 2 * a
      
      object@H[,idx] * (up/(down + eps))
    }, mc.cores = num_cores))
  } else {
    # Single core processing with optimized operations
    a <- WtW %*% object@H + H_pos + gamma_matrix
    b <- WtX
    c <- H_neg
    
    up <- b + sqrt(b * b + 4 * a * c)
    down <- 2 * a
    H_new <- object@H * (up/(down + eps))
  }
  
  # Return as dense matrix
  as(H_new, "dgeMatrix")
}

#' Update W and H Matrices
#'
#' @param object CSFNMF object
#' @param theta Convergence threshold
#' @param verbose Show progress
#' @param num_cores Number of cores for parallel processing
#' @return Updated object with optimization results
#' @importFrom progress progress_bar
#' @keywords internal
update_wh <- function(object, 
                      theta, 
                      verbose = TRUE, 
                      num_cores = 1) {
  
  # Initialize convergence monitoring
  loss <- numeric()
  loss_components <- list()
  loss[1] <- calculate_loss(object)
  loss_components[[1]] <- calculate_loss(object, divided = TRUE)
  
  # Initialize progress monitoring
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "Iteration :current/:total [:bar] :percent eta: :eta\n",
      total = object@parameters$max_iter
    )
  }
  
  # Initialize optimization state
  iter <- 1
  stop_cond <- 1
  best_state <- list(
    loss = loss[1],
    W = object@W,
    H = object@H
  )
  
  while (stop_cond > theta && 
         iter < object@parameters$max_iter && 
         loss[iter] > 0) {
    
    if (verbose) pb$tick()
    
    # Update matrices
    object@W <- update_w(object)
    object@H <- update_h(object, num_cores)
    
    # Calculate new loss
    loss[iter + 1] <- calculate_loss(object)
    
    # Update convergence condition
    stop_cond <- abs(loss[iter] - loss[iter + 1]) / abs(loss[iter])
    
    # Store best state if improved
    if (loss[iter + 1] < best_state$loss) {
      best_state <- list(
        loss = loss[iter + 1],
        W = object@W,
        H = object@H
      )
    }
    
    # Periodic garbage collection
    if (iter %% 10 == 0) gc()
    
    iter <- iter + 1
  }
  
  # Store optimization results
  object@results <- list(
    loss = loss[1:iter],
    loss_components = loss_components,
    iterations = iter - 1,
    convergence = stop_cond,
    final_loss = loss[iter],
    optimization_path = data.frame(
      iteration = 1:iter,
      loss = loss[1:iter],
      convergence = c(1, diff(loss[1:iter]))
    )
  )
  
  return(object)
}