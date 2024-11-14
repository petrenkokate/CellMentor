#' Train CSFNMF Model
#'
#' @param object CSFNMF object
#' @param pra_to_update Parameters to update ("alpha", "beta")
#' @param re_init Reinitialize matrices
#' @param num_update Number of updates
#' @param constSeed Random seed
#' @param verbose Show progress
#' @param theta Convergence threshold
#' @param num_cores Number of cores for parallel processing
#' @return Updated CSFNMF object
#' @export
TrainModelCSFNMF <- function(object,
                             pra_to_update = c("alpha", "beta"),
                             re_init = FALSE,
                             num_update = 10,
                             constSeed = 1,
                             verbose = TRUE,
                             theta = 1.0e-8,
                             num_cores = 1) {
  
  report <- create_reporter(verbose)
  train_object <- object@train_object
  
  # Initialize storage lists
  results <- list(
    loss = list(`no update` = train_object@results$loss),
    accuracy = list(`no update` = train_object@results$accuracy),
    WH = list(`no update` = list(W = train_object@W, H = train_object@H))
  )
  
  # Initial prediction
  singler_pred <- CSFnmfSingleR(train_object)
  report(sprintf("Initial accuracy: %.4f", train_object@results$accuracy))
  
  # Perform updates
  for (update in seq_len(num_update)) {
    report(sprintf("Update %d/%d", update, num_update))
    
    # Update hyperparameters
    new_paras <- with_memory_check({
      UpdateHyperpara(train_object, singler_pred, pra_to_update)
    })
    
    # Update helper matrices
    if ("alpha" %in% pra_to_update) {
      real_N <- train_object@constants@N
      train_object@constants@N <- new_paras$alpha * train_object@constants@N
      train_object@parameters$alpha <- sum(train_object@parameters$alpha)
    }
    
    if ("beta" %in% pra_to_update) {
      train_object@parameters$beta <- new_paras$beta
    }
    
    train_object@constants@Hconst <- calculate_const_for_h(train_object)
    
    # Reinitialize if requested
    if (re_init) {
      W0_H0 <- initialize_wh(train_object, train_object@parameters$init_method)
      train_object@W <- ensure_dgCMatrix(W0_H0$W)
      train_object@H <- ensure_dgCMatrix(W0_H0$H)
    }
    
    # Update W and H matrices
    train_object <- update_wh(train_object, theta, verbose)
    
    # Store results
    update_name <- sprintf("update_%d", update)
    results$WH[[update_name]] <- list(W = ensure_dgCMatrix(train_object@W),
                                      H = ensure_dgCMatrix(train_object@H))
    results$loss[[update_name]] <- train_object@results$loss
    
    # Restore alpha if updated
    if ("alpha" %in% pra_to_update) {
      train_object@parameters$alpha <- train_object@parameters$alpha * new_paras$alpha
      train_object@constants@N <- real_N
    }
    
    # Calculate accuracy
    h_project <- project_data(
      train_object@W,
      train_object@matrices@data,  
      seed = constSeed,
      num_cores = num_cores
    )
    
    accuracy_result <- calculate_accuracy(train_object, h_project, return_pred = TRUE)
    singler_pred <- accuracy_result$SingleRpred
    results$accuracy[[update_name]] <- accuracy_result$accuracy
    
    report(sprintf("Accuracy after update %d: %.4f", update, accuracy_result$accuracy))
  }
  
  # Select best update
  best_update <- names(which.max(unlist(results$accuracy)))
  report(sprintf("Best update: %s with accuracy %.4f", best_update, results$accuracy[[best_update]]))
  
  # Update object with best results
  train_object@results$accuracy <- results$accuracy[[best_update]]
  train_object@results$loss <- results$loss[[best_update]]
  best_W <- ensure_dgCMatrix(results$WH[[best_update]]$W)
  best_H <- ensure_dgCMatrix(results$WH[[best_update]]$H)
  
  train_object@W <- best_W
  train_object@H <- best_H
  object@W <- best_W
  
  # Calculate final H matrix
  h_combined <- cbind(object@H, h_project)
  object@H <- ensure_dgCMatrix(h_combined[, colnames(object@matrices@ref)])
  
  # Create update object with all results
  object@train_object <- methods::new(
    Class = "update_traincsfnmf",
    matrices = train_object@matrices,
    annotation = train_object@annotation,
    test_annotation = train_object@test_annotation,
    rank = train_object@parameters$rank,
    H = train_object@H,
    W = train_object@W,
    constants = train_object@constants,
    parameters = train_object@parameters,
    results = train_object@results,
    updates = list(
      WH_list = results$WH,
      loss_list = results$loss,
      accuracy_list = results$accuracy
    )
  )
  
  object
}

#' Update Hyperparameters
#'
#' @param train_object Training object
#' @param singler_pred SingleR predictions
#' @param pra_to_update Parameters to update
#' @return List of updated parameters
#' @keywords internal
UpdateHyperpara <- function(train_object, singler_pred, pra_to_update) {
  # Calculate correlation matrix
  correlation <- calculate_correlation_matrix(
    true_labels = train_object@test_annotation[singler_pred@rownames, "celltype"],
    pred_labels = singler_pred@listData[["labels"]]
  )
  
  # Update parameters
  list(
    beta = if ("beta" %in% pra_to_update) 
      update_beta(train_object, correlation) 
    else 
      train_object@parameters$beta,
    
    alpha = if ("alpha" %in% pra_to_update) 
      update_alpha(train_object, correlation) 
    else 
      train_object@parameters$alpha
  )
}

#' Calculate Correlation Matrix
#'
#' @param true_labels True cell type labels
#' @param pred_labels Predicted cell type labels
#' @return Correlation matrix
#' @keywords internal
calculate_correlation_matrix <- function(true_labels, pred_labels) {
  # Calculate initial correlation
  correlation <- table(true_labels, pred_labels)
  
  # Handle missing types
  missing_types <- setdiff(rownames(correlation), colnames(correlation))
  if (length(missing_types) > 0) {
    add_type <- matrix(
      0, 
      nrow = nrow(correlation), 
      ncol = length(missing_types),
      dimnames = list(rownames(correlation), missing_types)
    )
    correlation <- cbind(correlation, add_type)
  }
  
  # Zero out correct predictions
  to_zero <- outer(
    rownames(correlation),
    gsub('_[0-9]*', '', colnames(correlation)),
    "=="
  )
  correlation[to_zero] <- 0
  
  correlation
}

#' Update Beta Parameter
#'
#' @param train_object Training object
#' @param correlation Correlation matrix
#' @return Updated beta matrix
#' @keywords internal
update_beta <- function(train_object, correlation) {
  P <- train_object@constants@P
  num_pair_clusters <- ncol(P)
  diag_beta_new <- numeric(num_pair_clusters)
  beta_names <- character(num_pair_clusters)
  
  # Calculate new beta values
  for (col in seq_len(num_pair_clusters)) {
    cluster1 <- rownames(P)[which(P[,col] > 0)]
    cluster2 <- rownames(P)[which(P[,col] < 0)]
    parameter <- correlation[cluster1, cluster2] / max(sum(correlation[cluster1,]), 1)
    diag_beta_new[col] <- parameter
    beta_names[col] <- paste0(cluster1, "->", cluster2)
  }
  
  # Normalize beta values
  const_beta <- sum(train_object@parameters$beta) 
  diag_beta_new <- (diag_beta_new * const_beta) / sum(diag_beta_new)
  
  # Create beta matrix
  beta_new <- Matrix::Diagonal(x = diag_beta_new)
  colnames(beta_new) <- rownames(beta_new) <- beta_names
  
  beta_new
}

#' Update Alpha Parameter
#'
#' @param train_object Training object
#' @param correlation Correlation matrix
#' @return Updated alpha matrix
#' @importFrom matrixStats rowMedians
#' @keywords internal
update_alpha <- function(train_object, correlation) {
  num_cells <- ncol(train_object@H)
  num_cells_per_subtype <- table(train_object@annotation$celltype)
  diag_alpha_new <- numeric(num_cells)
  names(diag_alpha_new) <- colnames(train_object@H)
  
  # Process each cell type
  for (type in names(num_cells_per_subtype)) {
    # Get data for current type
    type_cells <- which(train_object@annotation$celltype == type)
    type_data <- train_object@H[, type_cells, drop = FALSE]
    
    # Calculate median and differences
    type_median <- matrixStats::rowMedians(as.matrix(type_data))
    matrix_median <- matrix(
      type_median,
      nrow = length(type_median),
      ncol = num_cells_per_subtype[type]
    )
    
    # Calculate distances
    diff_cells_median_type <- (type_data - matrix_median)^2
    type_SD_per_cell <- sqrt(colSums(diff_cells_median_type))
    
    # Calculate wrong predictions
    wrong_pred <- sum(correlation[type,]) / num_cells_per_subtype[type]
    
    # Set alpha values
    diag_alpha_new[colnames(type_data)] <- type_SD_per_cell * wrong_pred
    
  }
  
  # Normalize alpha values
  diag_alpha_new <- diag_alpha_new / sum(diag_alpha_new)
  
  # Create alpha matrix
  alpha_new <- Matrix::Diagonal(x = diag_alpha_new)
  colnames(alpha_new) <- rownames(alpha_new) <- names(diag_alpha_new)
  
  alpha_new
}