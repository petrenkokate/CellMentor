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
#' @importFrom tibble tibble
#' @keywords internal
initialize_cells_cluster <- function(data_matrix, k, labels){
  
  # browser()
  eps = .Machine$double.eps
  
  nFea = nrow(data_matrix)
  nSmp = ncol(data_matrix)
  
  W = matrix(0, nFea, k)
  H = matrix(0, k, nSmp)
  colnames(H) = colnames(data_matrix)
  rownames(W) = rownames(data_matrix)
  rownames(H)= colnames(W)= c(1:k)
  
  set.seed(1)
  
  num_clus = k
  data_for_clus = data_matrix
  data_for_clus = as.matrix(data_for_clus + eps) # KATE
  skmeans_clus = skmeans::skmeans(t(data_for_clus), num_clus, method="pclust")
  cell_clusters = skmeans_clus[["cluster"]]
  centroids = matrix(0, nrow= k, ncol=nrow(data_matrix))
  centroids[1:num_clus, ] = skmeans_clus[["prototypes"]]
  
  
  #For each cluster finding the main cell-type
  celltypeVSclus = table(labels, skmeans_clus[["cluster"]])
  mainType = tibble::tibble("skmeanClus" = c(1:k), 
                    "mainType" = rownames(celltypeVSclus)[apply(celltypeVSclus,2,which.max)])
  
  types = unique(labels)
  
  
  # Ensure that every cell type is represented by at least one cluster
  while (length(unique(mainType$mainType)) != length(types)){
    # Cell types that do not associated with any cluster
    non_main_types = types[which(!types %in% unique(mainType$mainType))]
    # Assign clusters to those cell types
    mainType$mainType[(num_clus - length(non_main_types) + 1):num_clus] = non_main_types
    num_clus = num_clus - length(non_main_types)
    
    for (clus in non_main_types){
      cells_name = names(which(labels == clus))
      clus_num = which(mainType$mainType == clus)
      cell_clusters[cells_name] = clus_num
      centroids[clus_num, ] = apply(data_matrix[, cells_name, drop = FALSE], 1, median)
    }
    
    data_for_clus = data_for_clus[, - which(labels[colnames(data_for_clus)] %in% non_main_types)]
    skmeans_clus = skmeans::skmeans(t(data_for_clus), num_clus, method="pclust")
    cell_clusters[names(skmeans_clus[["cluster"]])] = skmeans_clus[["cluster"]]
    celltypeVSclus = table(labels[names(skmeans_clus[["cluster"]])], skmeans_clus[["cluster"]])
    mainType$mainType[1:num_clus] = rownames(celltypeVSclus)[apply(celltypeVSclus,2,which.max)]
    centroids[1:num_clus, ] = skmeans_clus[["prototypes"]]
    
  }
  
  # Similarity <- function(cell, centroid) {lsa::cosine(cell, centroid)}
  Similarity <- function(cell, centroid) { #KATE
    if (all(cell == 0) || all(centroid == 0)) {
      return(0)  # Return zero if either vector is all zeros
    }
    lsa::cosine(cell, centroid)
  }
  
  # Each cell belonging to a cluster whose main cell-type does not match the cell's cell-type 
  # is associated to another cluster with the highest similarity that have the same main cell-type 
  missPredClus = colnames(data_matrix)[which(labels != mainType$mainType[cell_clusters])]
  newClus = as.integer(cell_clusters)
  names(newClus) = names(cell_clusters)
  
  relevantClus = similarMatrix = t(H)
  for (i in 1:k){
    similarMatrix[,i] = apply(data_matrix, 2, function(x){Similarity(x, centroids[i, ])})
  }
  
  # Only clusters that have the same main cell type
  for (i in rownames(relevantClus)){
    relevantClus[i, which(mainType$mainType == labels[i])] = 1
  }
  
  H = t(similarMatrix)
  
  similarMatrix = similarMatrix * relevantClus
  newClus[missPredClus] = unlist(apply(similarMatrix[missPredClus, ], 1, which.max))
  
  
  for (clus in names(table(newClus))){
    CellsName = names(which(newClus == clus))
    if (length(CellsName) == 1){
      W[, clus] = data_matrix[,CellsName]
    }
    else{
      W[, clus] = apply(data_matrix[,CellsName], 1, median)
    }
  }
  
  #Evoid zero elements
  W = pmax(W, eps)
  
  return(list("W" = Matrix(W, sparse=TRUE), "H" = Matrix(H, sparse=TRUE)))
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