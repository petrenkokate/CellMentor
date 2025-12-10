#' Create Training Object for CSFNMF
#'
#' @param ref_matrix Reference matrix
#' @param ref_celltype Reference cell types
#' @param data_matrix Data matrix
#' @param data_celltype Data cell types
#' @return TrainCSFNMF object
#' @keywords internal
#' @noRd
create_train_object <- function(ref_matrix,
                                ref_celltype,
                                data_matrix,
                                data_celltype) {
  # Create base object
  object <- methods::new("traincsfnmf")

  # Set matrices using new structure
  object@matrices <- methods::new(
    "RefDataList",
    ref = as(ref_matrix, "Matrix"),
    data = as(data_matrix, "Matrix")
  )

  # Set reference annotations
  names(ref_celltype) <- colnames(object@matrices@ref)
  object@annotation <- data.frame(
    celltype = ref_celltype,
    row.names = names(ref_celltype),
    stringsAsFactors = FALSE
  )

  # Set test annotations
  names(data_celltype) <- colnames(object@matrices@data)
  object@test_annotation <- data.frame(
    celltype = data_celltype,
    row.names = names(data_celltype),
    stringsAsFactors = FALSE
  )

  # Initialize other slots
  object@rank <- 0
  object@H <- Matrix(0)
  object@W <- Matrix(0)
  object@constants <- methods::new("helpmat")
  object@parameters <- list()
  object@results <- list()

  object
}

#' Divide Reference Data for Training
#'
#' @param object CSFNMF object
#' @param seed Random seed
#' @return Training object
#' @importFrom stats median
#' @importFrom Matrix summary rowSums
#' @importFrom data.table setDT
#' @keywords internal
#' @noRd
divide_reference_data <- function(object, seed = 1) {
  cell_id <- celltype <- .N <- . <- NULL

  # Convert to data.table and sample
  dt <- data.table::setDT(
    list(
      cell_id = rownames(object@annotation),
      celltype = object@annotation$celltype
    )
  )

  train_cells <- dt[, .(
    cell_id = sample(cell_id, ceiling(.N * 0.75))
  ), by = celltype][, cell_id]

  # Get data matrices
  data_matrix <- object@matrices@ref
  train_matrix <- data_matrix[, train_cells]

  # Handle zero genes using sparse matrix operations
  zero_genes <- which(Matrix::rowSums(train_matrix) == 0)
  if (length(zero_genes) > 0) {
    expr_matrix <- data_matrix[zero_genes, , drop = FALSE]
    nz <- Matrix::summary(expr_matrix)

    new_cells <- tapply(
      colnames(data_matrix)[nz$j],
      factor(nz$i, levels = seq_len(length(zero_genes))),
      function(x) if (length(x) > 0) sample(x, 1)
    )

    train_cells <- unique(c(train_cells, unlist(new_cells)))
  }

  # Create final matrices
  test_cells <- setdiff(colnames(data_matrix), train_cells)
  train_matrix <- data_matrix[, train_cells]
  test_matrix <- data_matrix[, test_cells]

  # Get cell types
  train_types <- object@annotation[train_cells, "celltype"]
  test_types <- object@annotation[test_cells, "celltype"]

  # Create training object
  create_train_object(
    train_matrix,
    train_types,
    test_matrix,
    test_types
  )
}

#' Project Data onto NMF Basis Matrix
#'
#' @description
#' Projects new data onto the learned basis matrix (W) using non-negative 
#' least squares (NNLS).
#' This function is used to obtain cell-type signatures (H matrix) for 
#' new query data
#' using the gene weights (W matrix) learned during training. The projection 
#' is performed
#' in chunks to manage memory efficiently, with optional parallel processing.
#'
#' @param W Basis matrix (genes × rank) containing learned gene weights
#' @param X Data matrix (genes × cells) to be projected. Must have same number 
#'          of genes (rows) as W
#' @param seed Random seed for reproducibility (default: 1)
#' @param num_cores Number of cores for parallel processing (default: 1).
#'                 If > 1, processing is parallelized across chunks
#' @param chunk_size Number of cells to process in each chunk (default: 1000).
#'                  Smaller chunks use less memory but may be slower
#' @param verbose Logical; whether to show progress bar (default: TRUE)
#'
#' @return A Matrix object (rank × cells) containing the projection 
#'        coefficients.
#'         The rows correspond to factors (rank) and columns to cells.
#'         Additional processing information is stored in attributes:
#'         - num_chunks: Number of chunks processed
#'         - chunk_size: Size of chunks used
#'         - num_cores: Number of cores used
#'
#' @details
#' The projection is performed using non-negative least squares (NNLS) to solve
#' the optimization problem: min ||X - WH||² subject to H >= 0, for each cell
#' in the input matrix X. The resulting H matrix contains the cell-type 
#' signatures for the query data.
#'
#' For memory efficiency, cells are processed in chunks. The chunk_size 
#' parameter can be adjusted based on available memory. Parallel processing can 
#'be enabled by setting num_cores > 1.
#'
#' @examples
#' # Minimal, fast example (no external data)
#' set.seed(1)
#'
#' # Dimensions
#' genes <- paste0("Gene", seq_len(50))
#' k     <- 3    # rank
#' cells <- 10
#'
#' # Non-negative basis W (genes x k)
#' W_ex <- matrix(abs(rnorm(length(genes) * k, sd = 0.5)),
#'                nrow = length(genes), ncol = k,
#'                dimnames = list(genes, paste0("k", seq_len(k))))
#'
#' # Generate a non-negative H_true and synthetic data X = W * H + noise
#' H_true <- matrix(abs(rnorm(k * cells, sd = 0.5)), nrow = k, ncol = cells)
#' X_ex   <- W_ex %*% H_true + matrix(rexp(length(genes) * cells, rate = 20),
#'                                    nrow = length(genes), ncol = cells,
#'                                    dimnames = list(genes, paste0("cell", seq_len(cells))))
#'
#' # Project (rank x cells)
#' H_est <- project_data(
#'   W = W_ex,
#'   X = X_ex,
#'   num_cores = 1,     # keep examples fast & deterministic
#'   chunk_size = 5,
#'   verbose = FALSE
#' )
#'
#' dim(H_est)           # should be k x cells
#' @importFrom nnls nnls
#' @importFrom BiocParallel bplapply SnowParam bpnworkers
#' @importFrom parallel detectCores
#' @importFrom progress progress_bar
#'
#' @export
project_data <- function(W, X, seed = 1, num_cores = 1,
                         chunk_size = 1000, verbose = TRUE) {
  k <- ncol(W)
  n_cells <- ncol(X)

  # Handle chunk_size
  if (is.null(chunk_size)) {
    chunk_size <- n_cells # one chunk for all cells
  } else if (chunk_size > n_cells) {
    chunk_size <- n_cells
    warning("chunk_size was larger than number of cells, using n_cells instead")
  }

  # Create progress bar if verbose
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "Projecting [:bar] :percent eta: :eta",
      total = ceiling(n_cells / chunk_size)
    )
  }

  # Function to process a chunk of cells
  process_chunk <- function(cell_indices) {
    # Pre-allocate matrix for chunk results
    chunk_size <- length(cell_indices)
    h_chunk <- matrix(0, nrow = k, ncol = chunk_size)

    # Process each cell in chunk
    for (i in seq_len(chunk_size)) {
      h_chunk[, i] <- nnls::nnls(W, X[, cell_indices[i]])$x
    }

    h_chunk
  }

  # Split cells into chunks
  cell_chunks <- split(
    seq_len(n_cells),
    ceiling(seq_len(n_cells) / chunk_size)
  )

  # Process chunks in parallel if requested
  if (num_cores > 1) {
    num_cores <- min(num_cores, parallel::detectCores(), length(cell_chunks))
    
    # Create BiocParallel backend
    BPPARAM <- BiocParallel::SnowParam(workers = num_cores)
    
    # Process chunks in parallel
    # Note: BiocParallel handles variable export automatically
    results <- BiocParallel::bplapply(
      cell_chunks, 
      function(chunk) {
        chunk_result <- process_chunk(chunk)
        if (verbose) pb$tick()
        chunk_result
      },
      BPPARAM = BPPARAM
    )
  } else {
    # Process sequentially
    results <- lapply(cell_chunks, function(chunk) {
      chunk_result <- process_chunk(chunk)
      if (verbose) pb$tick()
      chunk_result
    })
  }

  # Combine results
  h_proj <- do.call(cbind, results)

  # Set dimensions
  dimnames(h_proj) <- list(
    seq_len(k),
    colnames(X)
  )

  # Add attributes for debugging/monitoring
  attr(h_proj, "processing_info") <- list(
    num_chunks = length(cell_chunks),
    chunk_size = chunk_size,
    num_cores = num_cores
  )

  h_proj
}

#' Helper function for optimal chunk size determination
#'
#' @param n_cells Number of cells
#' @param available_memory Available memory in MB
#' @return Optimal chunk size
#' @keywords internal
#' @noRd
determine_chunk_size <- function(n_cells, available_memory = 1000) {
  # Estimate memory per cell (adjust based on your data)
  mem_per_cell <- 0.1 # MB

  # Calculate optimal chunk size
  chunk_size <- min(
    floor(available_memory / mem_per_cell),
    n_cells,
    1000 # Maximum chunk size
  )

  # Ensure minimum size
  max(chunk_size, 100)
}
