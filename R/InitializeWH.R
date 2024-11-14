#' Initialize W and H Matrices
#'
#' @description
#' Initializes W and H matrices using various methods
#'
#' @param object CSFNMF object
#' @param method Initialization method ("uniform", "regulated", "NNDSVD", "skmeanGenes", "skmeanCells")
#' @param seed Random seed
#' @return List containing W and H matrices
#' @importFrom Matrix Matrix
#' @importFrom sparsesvd sparsesvd
#' @importFrom cluster fanny
#' @importFrom skmeans skmeans
#' @export
initialize_wh <- function(object, method = "NNDSVD", seed = 1) {
  set.seed(seed)
  
  data_matrix <- as(object@matrices@ref, "Matrix")
  k <- object@parameters$rank  # Changed from object@rank@k
  celltype <- object@annotation$celltype
  
  # Removed celltype.code as it's not in new structure
  names(celltype) <- rownames(object@annotation)
  
  # Validate k
  max_possible_k <- min(nrow(data_matrix), ncol(data_matrix))
  if (k > max_possible_k) {
    warning(sprintf("Reducing k from %d to %d due to matrix dimensions", k, max_possible_k))
    k <- max_possible_k
  }
  
  switch(method,
         "uniform" = initialize_uniform(data_matrix, k),
         "regulated" = initialize_regulated(data_matrix, k, celltype),  # Changed from celltype.code
         "NNDSVD" = initialize_nndsvd(data_matrix, k),
         "skmeanGenes" = initialize_genes_cluster(data_matrix, k, celltype),
         "skmeanCells" = initialize_cells_cluster(data_matrix, k, celltype),
         stop("Unknown initialization method"))
}

#' Initialize Matrices Using Uniform Distribution
#'
#' @param data_matrix Data matrix
#' @param k Rank
#' @return List containing W and H matrices
#' @keywords internal
initialize_uniform <- function(data_matrix, k) {
  eps <- .Machine$double.eps
  nFea <- nrow(data_matrix)
  nSmp <- ncol(data_matrix)
  max_val <- max(data_matrix, na.rm = TRUE)
  
  # Initialize matrices using Matrix package
  H <- Matrix(
    data = runif(n = nSmp * k, min = eps, max = max_val),
    nrow = k,
    ncol = nSmp,
    sparse = TRUE
  )
  
  W <- Matrix(
    data = runif(n = nFea * k, min = eps, max = max_val),
    nrow = nFea,
    ncol = k,
    sparse = TRUE
  )
  
  # Set dimensions
  colnames(H) <- colnames(data_matrix)
  rownames(W) <- rownames(data_matrix)
  rownames(H) <- colnames(W) <- seq_len(k)
  
  list(W = W, H = H)
}

#' Initialize Matrices Using Up/Down Regulated Method
#'
#' @param data_matrix Data matrix
#' @param k Rank
#' @param labelcode Cell type codes
#' @return List containing W and H matrices
#' @keywords internal
initialize_regulated <- function(data_matrix, k, labelcode) {
  eps <- .Machine$double.eps
  nFea <- nrow(data_matrix)
  nSmp <- ncol(data_matrix)
  
  # Initialize H matrix
  H <- Matrix(
    data = runif(n = nSmp * k, min = eps),
    nrow = k,
    ncol = nSmp,
    sparse = TRUE
  )
  
  # Update H based on cell types
  types <- names(table(labelcode))
  for (i in seq_len(k)) {
    type_cells <- which(labelcode == types[i])
    if (length(type_cells) > 0) {
      H[i, type_cells] <- H[i, type_cells] + sum(H) / nSmp
    }
  }
  
  # Initialize W matrix
  W <- Matrix(
    data = runif(n = nFea * k, min = eps),
    nrow = nFea,
    ncol = k,
    sparse = TRUE
  )
  
  # Set dimensions
  colnames(H) <- colnames(data_matrix)
  rownames(W) <- rownames(data_matrix)
  rownames(H) <- colnames(W) <- seq_len(k)
  
  list(W = W, H = pmax(H, 0))
}

#' Initialize Using NNDSVD Method
#'
#' @param data_matrix Input matrix
#' @param k Rank
#' @return List containing W and H matrices
#' @importFrom Matrix sparseMatrix
#' @importFrom sparsesvd sparsesvd
#' @importFrom irlba irlba
#' @keywords internal
initialize_nndsvd <- function(data_matrix, k) {
  # Add safety checks
  if (k < 1) stop("k must be positive")
  
  # Get actual rank from SVD
  svd_result <- tryCatch({
    svd_dec <- sparsesvd::sparsesvd(data_matrix, rank = k)
    actual_k <- ncol(svd_dec$u)  # Get actual number of components returned
    if (actual_k < k) {
      warning(sprintf("sparsesvd returned only %d components instead of requested %d. irlba was used instead.", 
                      actual_k, k))
      # Fall back to irlba for full k components
      result <- irlba::irlba(data_matrix, nv = k)
      list(
        u = result$u,
        d = result$d,
        v = result$v
      )
    } else {
      svd_dec
    }
  }, error = function(e) {
    warning("sparsesvd failed, falling back to irlba")
    result <- irlba::irlba(data_matrix, nv = k)
    list(
      u = result$u,
      d = result$d,
      v = result$v
    )
  })
  
  U <- svd_result$u
  Sig <- svd_result$d
  V <- svd_result$v
  
  # Verify dimensions
  stopifnot(
    ncol(U) >= k,
    length(Sig) >= k,
    ncol(V) >= k
  )
  
  # Initialize matrices
  nFea <- nrow(data_matrix)
  nSmp <- ncol(data_matrix)
  W <- Matrix(0, nFea, k, sparse = TRUE)
  H <- Matrix(0, k, nSmp, sparse = TRUE)
  
  # First singular triplet
  W[,1] <- sqrt(Sig[1]) * abs(U[,1])
  H[1,] <- sqrt(Sig[1]) * abs(t(V[,1]))
  
  # Process remaining singular triplets
  for (i in 2:k) {
    u <- U[,i]
    v <- V[,i]
    
    # Calculate positive and negative parts
    upos <- pmax(u, 0)
    uneg <- abs(pmin(u, 0))
    vpos <- pmax(v, 0)
    vneg <- abs(pmin(v, 0))
    
    # Calculate norms with safety
    norms <- list(
      upos = sqrt(sum(upos^2)) + .Machine$double.eps,
      uneg = sqrt(sum(uneg^2)) + .Machine$double.eps,
      vpos = sqrt(sum(vpos^2)) + .Machine$double.eps,
      vneg = sqrt(sum(vneg^2)) + .Machine$double.eps
    )
    
    # Calculate terms
    termsq_pos <- norms$upos * norms$vpos
    termsq_neg <- norms$uneg * norms$vneg
    
    if (termsq_pos >= termsq_neg) {
      W[,i] <- sqrt(Sig[i] * termsq_pos) * upos / norms$upos
      H[i,] <- sqrt(Sig[i] * termsq_pos) * vpos / norms$vpos
    } else {
      W[,i] <- sqrt(Sig[i] * termsq_neg) * uneg / norms$uneg
      H[i,] <- sqrt(Sig[i] * termsq_neg) * vneg / norms$vneg
    }
  }
  
  # Handle small values
  eps <- .Machine$double.eps
  W[W < eps] <- eps
  H[H < eps] <- eps
  
  # Set dimensions
  colnames(H) <- colnames(data_matrix)
  rownames(W) <- rownames(data_matrix)
  rownames(H) <- colnames(W) <- seq_len(k)
  
  list(
    W = ensure_dgCMatrix(W),
    H = ensure_dgCMatrix(H)
  )
}

#' Initialize Matrices Using Gene Clustering
#'
#' @param data_matrix Data matrix
#' @param k Rank
#' @param labels Cell type labels
#' @return List containing W and H matrices
#' @importFrom cluster fanny
#' @keywords internal
initialize_genes_cluster <- function(data_matrix, k, labels) {
  # Perform fuzzy clustering
  clusters <- cluster::fanny(data_matrix, k = k)
  W <- clusters$membership
  H <- t(W) %*% data_matrix
  
  # Process by cell type
  cell_types <- unique(labels)
  for (cell_type in cell_types) {
    cells <- colnames(data_matrix)[labels == cell_type]
    if (length(cells) > 0) {
      H[, cells] <- apply(H[, cells, drop = FALSE], 1, median)
    }
  }
  
  # Set dimensions
  colnames(H) <- colnames(data_matrix)
  rownames(W) <- rownames(data_matrix)
  rownames(H) <- colnames(W) <- seq_len(k)
  
  list(W = as(W, "Matrix"), H = as(H, "Matrix"))
}

#' Initialize Matrices Using Cell Clustering
#'
#' @param data_matrix Data matrix
#' @param k Rank
#' @param labels Cell type labels
#' @return List containing W and H matrices
#' @importFrom skmeans skmeans
#' @keywords internal
initialize_cells_cluster <- function(data_matrix, k, labels) {
  eps <- .Machine$double.eps
  
  # Initialize matrices
  nFea <- nrow(data_matrix)
  nSmp <- ncol(data_matrix)
  W <- matrix(0, nFea, k)
  H <- matrix(0, k, nSmp)
  
  # Perform spherical k-means clustering
  clusters <- skmeans::skmeans(
    t(as.matrix(data_matrix + eps)),
    k,
    method = "pclust"
  )
  
  # Get cell clusters and centroids
  cell_clusters <- clusters$cluster
  centroids <- matrix(0, nrow = k, ncol = nFea)
  centroids[1:k,] <- clusters$prototypes
  
  # Process cell type relationships
  type_cluster_table <- table(labels, cell_clusters)
  main_types <- tibble::tibble(
    skmean_clus = seq_len(k),
    main_type = rownames(type_cluster_table)[apply(type_cluster_table, 2, which.max)]
  )
  
  # Ensure all cell types are represented
  types <- unique(labels)
  if (length(unique(main_types$main_type)) != length(types)) {
    # Handle missing cell types
    missing_types <- setdiff(types, unique(main_types$main_type))
    remaining_clusters <- k - length(missing_types)
    
    for (type in missing_types) {
      cells <- colnames(data_matrix)[labels == type]
      cluster_num <- which(main_types$main_type == main_types$main_type[remaining_clusters])
      cell_clusters[cells] <- cluster_num
      centroids[cluster_num,] <- apply(data_matrix[,cells], 1, median)
      remaining_clusters <- remaining_clusters - 1
    }
  }
  
  # Calculate similarity matrix
  H <- calculate_similarity_matrix(data_matrix, centroids, cell_clusters)
  
  # Calculate W matrix
  for (cluster in seq_len(k)) {
    cells <- names(cell_clusters)[cell_clusters == cluster]
    if (length(cells) == 1) {
      W[, cluster] <- data_matrix[,cells]
    } else {
      W[, cluster] <- apply(data_matrix[,cells], 1, median)
    }
  }
  
  # Clean up near-zero elements
  W <- pmax(W, eps)
  
  # Set dimensions
  colnames(H) <- colnames(data_matrix)
  rownames(W) <- rownames(data_matrix)
  rownames(H) <- colnames(W) <- seq_len(k)
  
  list(W = as(W, "Matrix"), H = as(H, "Matrix"))
}

#' Calculate Similarity Matrix
#'
#' @param data_matrix Data matrix
#' @param centroids Centroid matrix
#' @param cell_clusters Cell cluster assignments
#' @return Similarity matrix
#' @importFrom lsa cosine
#' @keywords internal
calculate_similarity_matrix <- function(data_matrix, centroids, cell_clusters) {
  k <- nrow(centroids)
  n_cells <- ncol(data_matrix)
  
  # Initialize similarity matrix
  H <- matrix(0, k, n_cells)
  
  # Calculate similarities
  for (i in seq_len(k)) {
    H[i,] <- vapply(seq_len(n_cells), function(j) {
      cell <- data_matrix[,j]
      if (all(cell == 0) || all(centroids[i,] == 0)) {
        return(0)
      }
      lsa::cosine(cell, centroids[i,])
    }, numeric(1))
  }
  
  H
}