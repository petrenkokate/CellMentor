#' Calculate Loss Function
#'
#' @param object CSFNMF object
#' @param divided Logical: return components separately
#' @return Loss value or list of components
#' @importFrom Matrix norm
#' @keywords internal
calculate_loss <- function(object, divided = FALSE) {
  X <- object@matrices@ref
  W <- object@W
  H <- object@H
  const.list <- object@parameters
  num.clusters <- length(table(object@annotation$celltype))
  
  # Calculate components with memory optimization
  error <- with_memory_check({
    norm(as(X - W %*% H, "matrix"), type = "F")^2
  })
  
  within_var <- with_memory_check({
    H - H %*% object@constants@A
  })
  
  sw <- with_memory_check({
    sum(diag(as.matrix(within_var) %*% as.matrix(object@constants@N) %*% t(as.matrix(within_var))))
  })
  
  between_var <- with_memory_check({
    H %*% object@constants@BP
  })
  
  sb <- with_memory_check({
    sum(diag(as.matrix(between_var) %*% as.matrix(const.list[["beta"]]) %*% t(as.matrix(between_var))))
  })
  
  reg.h <- sum(H)
  reg.w <- sum(diag(t(as.matrix(W)) %*% as.matrix(W) %*% as.matrix(object@constants@M)))
  
  if (divided) {
    list(
      error = 0.5 * error,
      sw = (const.list[["alpha"]]/(2*num.clusters)) * sw,
      sb = (1/(2*(num.clusters-1)*num.clusters)) * sb,
      reg.h = const.list[["gamma"]] * reg.h,
      reg.w = (const.list[["delta"]]/2) * reg.w
    )
  } else {
    0.5 * error + 
      (const.list[["alpha"]]/(2*num.clusters)) * sw - 
      (1/(2*(num.clusters-1)*num.clusters)) * sb + 
      const.list[["gamma"]] * reg.h + 
      (const.list[["delta"]]/2) * reg.w
  }
}

#' Update W Matrix
#'
#' @param object CSFNMF object
#' @return Updated W matrix
#' @importFrom Matrix Diagonal
#' @keywords internal
update_w <- function(object) {
  eps <- .Machine$double.eps
  
  # Calculate components in chunks to optimize memory
  up <- with_memory_check({
    object@matrices@ref %*% t(object@H)
  })
  
  down <- with_memory_check({
    object@W %*% (
      object@H %*% t(object@H) + 
        object@parameters[["delta"]] * object@constants@M
    )
  })
  
  # Update W
  W_new <- object@W * (up/(down + eps))
  
  as(W_new, "Matrix")
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
  
  # Set up parallel processing
  if (num_cores > 1) {
    num_cores <- min(num_cores, parallel::detectCores())
  }
  
  # Calculate components
  parts <- with_memory_check({
    list(
      a = crossprod(object@W, object@W) %*% object@H + 
        object@H %*% object@constants@Hconst[["positive"]] + 
        matrix(object@parameters[["gamma"]], nrow = nrow(object@H), 
               ncol = ncol(object@H)),
      b = crossprod(object@W, object@matrices@ref),
      c = object@H %*% object@constants@Hconst[["negative"]]
    )
  })
  
  # Update H in parallel if requested
  if (num_cores > 1) {
    # Split matrix into chunks
    chunks <- split(seq_len(ncol(object@H)), 
                    cut(seq_len(ncol(object@H)), num_cores))
    
    # Process chunks in parallel
    H_new <- do.call(cbind, parallel::mclapply(chunks, function(idx) {
      up <- parts$b[,idx] + sqrt(parts$b[,idx] * parts$b[,idx] + 
                                   4 * parts$a[,idx] * parts$c[,idx])
      down <- 2 * parts$a[,idx]
      object@H[,idx] * (up/(down + eps))
    }, mc.cores = num_cores))
  } else {
    # Single core processing
    up <- parts$b + sqrt(parts$b * parts$b + 4 * parts$a * parts$c)
    down <- 2 * parts$a
    H_new <- object@H * (up/(down + eps))
  }
  
  as(H_new, "Matrix")
}

#' Update W and H Matrices
#'
#' @param object CSFNMF object
#' @param theta Convergence threshold
#' @param verbose Show progress
#' @param num_cores Number of cores for parallel processing
#' @return Updated object
#' @importFrom progress progress_bar
#' @keywords internal
update_wh <- function(object, theta, verbose = TRUE, num_cores = 1) {
  # Initialize progress bar
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "Iteration :current/:total [:bar] :percent eta: :eta",
      total = object@parameters$max_iter
    )
  }
  
  # Initialize variables
  loss <- numeric()
  loss[1] <- calculate_loss(object)
  iter <- 1
  stop_cond <- 1
  
  while (stop_cond > theta && iter < object@parameters$max_iter && loss[iter] > 0) {
    if (verbose) pb$tick()
    
    # Update matrices with memory optimization
    with_memory_check({
      # Update W
      object@W <- update_w(object)
      
      # Update H
      object@H <- update_h(object, num_cores)
      
      # Calculate new loss
      loss[iter + 1] <- calculate_loss(object)
      
      # Check convergence
      stop_cond <- loss[iter] - loss[iter + 1]
      iter <- iter + 1
    })
    
    # Garbage collection every few iterations
    if (iter %% 5 == 0) gc()
  }
  
  # Store loss in results list instead of directly in object
  object@results$loss <- loss
  object
}