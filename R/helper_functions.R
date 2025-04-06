#' Calculate Standard Beta Matrix
#'
#' @description
#' Creates standardized beta matrix for NMF calculations
#'
#' @param object CSFNMF object
#' @param const Constant value for diagonal (default = 1)
#' @return Beta matrix
#' @importFrom Matrix Diagonal
#' @keywords internal
standard_beta <- function(object, const = 1) {
  # Calculate dimensions
  num_clusters <- length(table(object@annotation$celltype))
  dim_beta <- num_clusters * (num_clusters - 1)
  
  # Create diagonal matrix
  beta <- Matrix::Diagonal(
    n = dim_beta,
    x = const/dim_beta
  )
  
  # Convert to sparse matrix
  as(beta, "Matrix")
}

#' Calculate N Matrix
#'
#' @param celltype Cell type vector
#' @return N matrix
#' @importFrom Matrix Diagonal
#' @keywords internal
calculate_n_matrix <- function(celltype) {
  # Calculate inverse counts
  type_counts <- table(celltype)
  cell_weights <- 1/type_counts[celltype]
  
  # Create diagonal matrix with proper names
  N <- Matrix::Diagonal(x = as.numeric(cell_weights))
  
  # Set row and column names
  colnames(N) <- rownames(N) <- names(celltype)
  
  as(N, "Matrix") 
}

#' Calculate Help Matrices
#'
#' @description
#' Calculates helper matrices for NMF optimization
#'
#' @param object CSFNMF object
#' @return Updated object with calculated matrices
#' @importFrom Matrix bdiag
#' @keywords internal
calculate_help_matrices <- function(object) {
  # Get cell type counts and info
  celltype <- object@annotation$celltype
  names(celltype) <- rownames(object@annotation)  # Make sure celltype has names
  counts <- table(celltype)
  num_clusters <- length(counts)
  k <- object@parameters$rank 
  
  # Calculate matrices

  # Calculate A matrix
  A <- calculate_a_matrix(counts)
  # Ensure A has names
  colnames(A) <- rownames(A) <- names(celltype)
  
  # Calculate B matrix
  B <- calculate_b_matrix(counts)
  # Ensure B has proper row names
  rownames(B) <- names(celltype)
  colnames(B) <- names(counts)
  
  # Calculate P matrix
  P <- calculate_p_matrix(num_clusters)
  rownames(P) <- names(counts)
  pair_names <- character(ncol(P))
  col_idx <- 1
  for (i in seq_len(num_clusters)) {
    for (j in setdiff(seq_len(num_clusters), i)) {
      pair_names[col_idx] <- sprintf("%s-%s", names(counts)[i], names(counts)[j])
      col_idx <- col_idx + 1
    }
  }
  colnames(P) <- pair_names
  
  # Calculate M matrix
  M <- as(matrix(1, ncol = k, nrow = k) - 
            diag(x = 1, nrow = k, ncol = k), 
          "Matrix")
  rownames(M) <- colnames(M) <- paste0("factor_", seq_len(k))
  
  # Calculate N matrix with names
  N <- calculate_n_matrix(celltype)
  
  # Calculate BP and BP_posneg
  BP <- B %*% P
  rownames(BP) <- rownames(B)
  colnames(BP) <- colnames(P)
  
  BP_posneg <- calculate_bp_posneg(B, P)
  
  # Create help matrices object
  const_matrices <- methods::new(
    Class = "helpmat",
    A = A,
    B = B,
    P = P,
    M = M,
    N = N,
    BP = BP,
    BP_posneg = BP_posneg,
    Hconst = list()
  )
  
  const_matrices
}

#' Calculate A Matrix
#'
#' @param counts Cell type counts
#' @return A matrix
#' @keywords internal
calculate_a_matrix <- function(counts) {
  # Create block diagonal matrix
  A <- Matrix::bdiag(
    lapply(counts, function(x) {
      matrix(1/x, x, x)
    })
  )
  
  as(A, "Matrix")
}

#' Calculate B Matrix
#'
#' @param counts Cell type counts
#' @return B matrix
#' @keywords internal
calculate_b_matrix <- function(counts) {
  num_cell <- sum(counts)
  num_clusters <- length(counts)
  
  # Create sparse matrix
  B <- Matrix(0, num_cell, num_clusters)
  counter <- 0
  
  for (i in seq_along(counts)) {
    num_type <- counts[[i]]
    range <- c(counter + 1, counter + num_type)
    B[range[1]:range[2], i] <- 1/num_type
    counter <- counter + num_type
  }
  
  as(B, "Matrix")
}

#' Calculate P Matrix
#'
#' @param num_clusters Number of clusters
#' @return P matrix
#' @keywords internal
calculate_p_matrix <- function(num_clusters) {
  if (num_clusters == 1) {
    return(1)
  }
  
  if (num_clusters == 2) {
    P <- matrix(1, nrow = 2, ncol = 2)
    P[2,1] <- P[1,2] <- -1
    return(as(P, "Matrix"))
  }
  
  # Create P matrix for more than 2 clusters
  total_cols <- num_clusters * (num_clusters - 1)
  P <- Matrix(0, num_clusters, total_cols)
  
  # Fill P matrix
  col_idx <- 1
  for (c in seq_len(num_clusters)) {
    other_clusters <- setdiff(seq_len(num_clusters), c)
    for (o in other_clusters) {
      P[c, col_idx] <- 1
      P[o, col_idx] <- -1
      col_idx <- col_idx + 1
    }
  }
  
  as(P, "Matrix")
}

#' Calculate BP Positive/Negative Components
#'
#' @param B B matrix
#' @param P P matrix
#' @return List of positive and negative components
#' @keywords internal
calculate_bp_posneg <- function(B, P) {
  # Split P into positive and negative components
  P_pos <- P
  P_pos@x[P_pos@x < 0] <- 0
  
  P_neg <- P
  P_neg@x[P_neg@x > 0] <- 0
  P_neg@x <- abs(P_neg@x)
  
  # Calculate components
  list(
    positive = B %*% P_pos,
    negative = B %*% P_neg
  )
}

#' Calculate Constants for H
#'
#' @param object CSFNMF object
#' @param beta Optional beta matrix
#' @return List of constants for H updates
#' @keywords internal
calculate_const_for_h <- function(object, beta = NULL) {
  # Use provided beta or get from parameters
  if (is.null(beta)) {
    beta <- object@parameters$beta  # Changed from hyper_para
  }
  
  # Get components
  pos <- object@constants@BP_posneg$positive
  neg <- object@constants@BP_posneg$negative
  
  # Calculate matrix products
  P_pos_pos <- pos %*% beta %*% t(pos)
  P_pos_neg <- pos %*% beta %*% t(neg)
  P_neg_pos <- neg %*% beta %*% t(pos)
  P_neg_neg <- neg %*% beta %*% t(neg)
  
  # Calculate constants
  num_clusters <- length(table(object@annotation$celltype))
  sw_const <- object@parameters$alpha / num_clusters  # Changed from hyper_para
  sb_const <- 1/((num_clusters - 1) * num_clusters)
  
  A <- object@constants@A
  N <- object@constants@N
  
  # Calculate negative and positive constants
  neg_const <- sw_const * (N %*% A + A %*% N) + 
    sb_const * (P_pos_pos + P_neg_neg)
  
  pos_const <- sw_const * (N + A %*% N %*% A) + 
    sb_const * (P_neg_pos + P_pos_neg)
  
  list(
    positive = pos_const,
    negative = neg_const
  )
}

#' Calculate Accuracy
#'
#' @param train_object Training object
#' @param h_project Projected H matrix
#' @param seed Random seed
#' @param return_pred Return prediction details
#' @return Accuracy or list with accuracy and predictions
#' @importFrom MLmetrics Accuracy
#' @keywords internal
calculate_accuracy <- function(train_object, h_project, 
                               seed = 1, return_pred = FALSE) {
  # Get predictions
  singler_pred <- CSFnmfSingleR(train_object, h_project)
  
  # Get labels
  pred_labels <- singler_pred@listData[["labels"]]
  true_labels <- train_object@test_annotation[singler_pred@rownames, "celltype"]
  
  # Clean labels
  pred_clean <- gsub('[0-9]', '', pred_labels)
  true_clean <- gsub('[0-9]', '', true_labels)
  
  # Calculate accuracy
  accuracy <- MLmetrics::Accuracy(pred_clean, true_clean)
  
  if (return_pred) {
    list(
      accuracy = accuracy,
      SingleRpred = singler_pred
    )
  } else {
    accuracy
  }
}

#' SingleR Prediction on CSFNMF
#'
#' @param train_object Training object
#' @param H_data Optional H matrix
#' @param de.n Number of genes for DE
#' @return SingleR predictions
#' @importFrom SingleR SingleR
#' @keywords internal
CSFnmfSingleR <- function(train_object, H_data = NULL, de.n = 50) {
  # Get reference projections
  H_ref <- as.data.frame(project_data(
    train_object@W,
    train_object@matrices@ref  # Changed from count.matrices
  ))
  
  k <- train_object@parameters$rank  # Changed from rank@k
  rownames(H_ref) <- seq_len(k)
  
  # Get data projections if not provided
  if (is.null(H_data)) {
    H_data <- as.data.frame(project_data(
      train_object@W,
      train_object@matrices@data  # Changed from count.matrices
    ))
  }
  
  # Run SingleR
  SingleR::SingleR(
    test = H_data,
    ref = H_ref[, rownames(train_object@annotation)],
    labels = train_object@annotation$celltype,
    de.n = de.n
  )
}

#' Calculate Normalized Mutual Information
#'
#' @param true_labels Vector of true labels
#' @param pred_labels Vector of predicted labels
#' @return NMI score between 0 and 1
#' @importFrom entropy entropy.empirical mi.empirical
#' @keywords internal
calculate_nmi <- function(true_labels, pred_labels) {
  
  NMI(true_labels, pred_labels)
  
}

#' Perform clustering on NMF factors
#'
#' @param h_matrix H matrix from NMF
#' @param n_clusters Number of clusters to find
#' @param method Clustering method ('kmeans', 'hclust', or 'leiden')
#' @return Vector of cluster assignments
#' @importFrom stats kmeans hclust cutree dist
#' @keywords internal
cluster_nmf <- function(h_matrix, n_clusters, method = "kmeans", resolution = 0.8) {
  # Transpose to have cells as rows
  data <- t(as.matrix(h_matrix))
  
  clusters <- switch(
    method,
    "kmeans" = {
      set.seed(1)  # for reproducibility
      km <- kmeans(data, centers = n_clusters)
      km$cluster
    },
    "hclust" = {
      hc <- hclust(dist(data))
      cutree(hc, k = n_clusters)
    },
    "snn" = {
      # Create Seurat object
      seur_obj <- CreateSeuratObject(counts = t(data))
      assay_v3 <- CreateAssayObject(
        counts = seur_obj[["RNA"]]$counts
      )
      
      seur_obj[["RNA"]] <- assay_v3
      seur_obj <- FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
      seur_obj <- FindNeighbors(seur_obj, dims = NULL, annoy.metric="cosine", verbose = FALSE)
      seur_obj <- FindClusters(seur_obj, resolution = resolution, algorithm = 1, verbose = FALSE)
      as.numeric(Idents(seur_obj))
    },
    stop("Unknown clustering method")
  )
  
  clusters
}

#' Calculate Performance Metrics using Clustering
#'
#' @param train_object Training object
#' @param h_project Projected H matrix
#' @param seed Random seed
#' @param return_pred Return prediction details
#' @param method Clustering method to use
#' @return List containing NMI score and optionally cluster assignments
#' @importFrom aricode NMI
#' @keywords internal
calculate_performance <- function(train_object, h_project, 
                                  seed = 1, return_pred = FALSE,
                                  method = "snn") {
  # Get number of true clusters
  n_clusters <- length(unique(train_object@test_annotation$celltype))
  
  # Perform clustering
  clusters <- cluster_nmf(h_project, n_clusters, method)
  
  # Get cell names in correct order
  cell_names <- colnames(train_object@matrices@data)
  
  # Get true labels
  true_labels <- train_object@test_annotation[cell_names, "celltype"]
  
  # Calculate NMI
  nmi_score <- aricode::NMI(true_labels, clusters)
  
  if (return_pred) {
    list(
      nmi = nmi_score,
      clusters = clusters,
      true_labels = true_labels
    )
  } else {
    nmi_score
  }
}

#' Update Training Object Performance with Clustering Results
#'
#' @param train_object Training object
#' @param h_project Projected H matrix 
#' @param num_cores Number of cores for parallel processing
#' @param method Clustering method to use
#' @return Updated training object with performance metrics
#' @keywords internal
update_training_performance <- function(train_object, h_project, 
                                        num_cores = 1, method = "kmeans") {
  # Calculate performance with clustering
  perf_result <- calculate_performance(
    train_object, 
    h_project, 
    return_pred = TRUE,
    method = method
  )
  
  # Update training object results
  train_object@results$nmi <- perf_result$nmi
  train_object@results$clusters <- perf_result$clusters
  train_object@results$true_labels <- perf_result$true_labels
  
  train_object
}