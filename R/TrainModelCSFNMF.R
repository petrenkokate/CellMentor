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
  full_object <- methods::new("traincsfnmf")
  full_object@matrices <- object@matrices
  full_object@annotation <- object@annotation
  
  full_object@parameters <- train_object@parameters
  full_object@constants <- train_object@constants
  full_object@W <- object@W
  full_object@H <- object@H
  # Initialize storage lists
  results <- list(
    loss = list(`no update` = train_object@results$loss),
    nmi = list(`no update` = train_object@results$nmi),
    WH = list(`no update` = list(W = full_object@W, H = full_object@H))
  )
  
  # Initial prediction
  h_project <- project_data(
    W = full_object@W,
    X = train_object@matrices@data,
    seed = constSeed,
    num_cores = num_cores
  )
  init_clusters <- calculate_performance(train_object, h_project, return_pred = T)
  report(sprintf("Initial NMI: %.4f", train_object@results$nmi))
  clusters <- init_clusters$clusters 
  # Perform updates
  for (update in seq_len(num_update)) {
    report(sprintf("Update %d/%d", update, num_update))
    
    # Update hyperparameters
    new_paras <- UpdateHyperpara(train_object, full_object, clusters, pra_to_update)

    
    # Update helper matrices
    if ("alpha" %in% pra_to_update) {
      real_N <- full_object@constants@N
      full_object@constants@N <- new_paras$alpha * full_object@constants@N
      full_object@parameters$alpha <- sum(full_object@parameters$alpha)
    }
    
    if ("beta" %in% pra_to_update) {
      full_object@parameters$beta <- new_paras$beta
    }
    
    full_object@constants@Hconst <- calculate_const_for_h(full_object)
    
    # Reinitialize if requested
    if (re_init) {
      W0_H0 <- initialize_wh(full_object, full_object@parameters$init_method)
      full_object@W <- ensure_dgCMatrix(W0_H0$W)
      full_object@H <- ensure_dgCMatrix(W0_H0$H)
    }
    
    # Update W and H matrices
    full_object <- update_wh(full_object, theta, verbose)
    
    # Store results
    update_name <- sprintf("update_%d", update)
    results$WH[[update_name]] <- list(W = ensure_dgCMatrix(full_object@W),
                                      H = ensure_dgCMatrix(full_object@H))
    results$loss[[update_name]] <- full_object@results$loss
    
    # Restore alpha if updated
    if ("alpha" %in% pra_to_update) {
      full_object@parameters$alpha <- full_object@parameters$alpha * new_paras$alpha
      full_object@constants@N <- real_N
    }
    
    # Calculate NMI
    h_project <- project_data(
      full_object@W,
      train_object@matrices@data,  
      seed = constSeed,
      num_cores = num_cores
    )
    
    nmi_result <- calculate_performance(train_object, h_project, return_pred = TRUE)
    clusters <- nmi_result$clusters
    results$nmi[[update_name]] <- nmi_result$nmi
    
    report(sprintf("NMI after update %d: %.4f", update, nmi_result$nmi))
  }
  
  # Select best update
  best_update <- names(which.max(unlist(results$nmi)))
  report(sprintf("Best update: %s with NMI %.4f", best_update, results$nmi[[best_update]]))
  
  # Update object with best results
  train_object@results$nmi <- results$nmi[[best_update]]
  train_object@results$loss <- results$loss[[best_update]]
  best_W <- ensure_dgCMatrix(results$WH[[best_update]]$W)
  best_H <- ensure_dgCMatrix(results$WH[[best_update]]$H)
  
  train_object@W <- best_W
  train_object@H <- best_H[, colnames(train_object@matrices@ref)]
  object@W <- best_W
  
  # Calculate final H matrix
  h_combined <- cbind(object@H, h_project)
  object@H <- best_H
  
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
      nmi_list = results$nmi
    )
  )
  
  object
}

#' Update Hyperparameters
#'
#' @param train_object Training object
#' @param clusters Clustering assignments
#' @param pra_to_update Parameters to update
#' @return List of updated parameters
#' @keywords internal
UpdateHyperpara <- function(train_object, full_object, clusters, pra_to_update) {
  # Calculate correlation matrix
  correlation <- calculate_correlation_matrix(
    true_labels = train_object@test_annotation$celltype,
    clusters = clusters
  )
  
  # Update parameters
  list(
    beta = if ("beta" %in% pra_to_update) 
      update_beta(full_object, correlation) 
    else 
      full_object@parameters$beta,
    
    alpha = if ("alpha" %in% pra_to_update) 
      update_alpha(full_object, correlation) 
    else 
      full_object@parameters$alpha
  )
}

#' Calculate Correlation Matrix
#'
#' @param true_labels True cell type labels
#' @param clusters Cluster assignments
#' @return Correlation matrix
#' @keywords internal
calculate_correlation_matrix <- function(true_labels, clusters) {
  # Calculate initial correlation
  correlation <- table(true_labels, clusters)
  
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
  
  # Normalize by row sums (optional: makes it proportional)
  correlation <- correlation / rowSums(correlation)
  
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
  
  # Calculate new beta values based on cluster correlations
  for (col in seq_len(num_pair_clusters)) {
    cluster1 <- rownames(P)[which(P[,col] > 0)]
    cluster2 <- rownames(P)[which(P[,col] < 0)]
    
    # Calculate separation between clusters
    separation <- mean(correlation[cluster1,]) - correlation[cluster1, cluster2]
    diag_beta_new[col] <- max(separation, 0)
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
  # Get all cell names from N matrix
  all_cells <- colnames(train_object@constants@N)
  
  # Create full-sized alpha matrix
  diag_alpha_new <- numeric(length(all_cells))
  names(diag_alpha_new) <- all_cells
  
  # Process each cell type
  for (type in names(table(train_object@annotation$celltype))) {
    type_cells <- which(train_object@annotation$celltype == type)
    
    # Your existing calculations for cells of this type
    type_data <- train_object@H[, type_cells, drop = FALSE]
    type_median <- matrixStats::rowMedians(as.matrix(type_data))
    matrix_median <- matrix(
      type_median,
      nrow = length(type_median),
      ncol = length(type_cells)
    )
    
    diff_cells_median_type <- (type_data - matrix_median)^2
    type_SD_per_cell <- sqrt(colSums(diff_cells_median_type))
    
    # Get wrong cluster proportion for this type
    wrong_cluster <- sum(correlation[type, ]) / length(type_cells)
    
    # Set values in full alpha vector
    diag_alpha_new[type_cells] <- type_SD_per_cell * wrong_cluster
  }
  
  # Normalize
  diag_alpha_new <- diag_alpha_new / sum(diag_alpha_new)
  
  # Create matrix
  alpha_new <- Matrix::Diagonal(x = diag_alpha_new)
  colnames(alpha_new) <- rownames(alpha_new) <- all_cells
  
  alpha_new
}