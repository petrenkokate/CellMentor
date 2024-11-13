#' Null coalescing operator
#'
#' @param x First value
#' @param y Default value
#' @return x if not NULL, otherwise y
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Memory-Safe Operation
#'
#' @description
#' Executes operation with memory monitoring and cleanup
#'
#' @param expr Expression to evaluate
#' @param threshold Memory threshold in MB
#' @return Result of expression
#' @importFrom utils memory.size
#' @keywords internal
with_memory_check <- function(expr, threshold = 10000) {
  # Initial cleanup
  gc()
  initial_mem <- memory_usage()
  
  # Evaluate expression
  result <- tryCatch({
    expr
  }, error = function(e) {
    gc()
    stop("Operation failed: ", e$message)
  })
  
  # Check memory usage
  final_mem <- memory_usage()
  mem_diff <- final_mem - initial_mem
  
  if (mem_diff > threshold) {
    warning(sprintf("Large memory increase: %.2f MB", mem_diff))
    gc()
  }
  
  result
}

#' Create Progress Manager
#'
#' @param total Total steps
#' @param title Progress bar title
#' @return Progress manager object
#' @importFrom progress progress_bar
#' @keywords internal
create_progress <- function(total, title = "Progress") {
  list(
    bar = progress::progress_bar$new(
      format = paste0(
        title, " [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
        "\n:message"
      ),
      total = total,
      width = 80
    ),
    start_time = Sys.time(),
    memory_start = memory_usage()
  )
}

#' Update Progress
#'
#' @param progress Progress manager
#' @param message Progress message
#' @param show_memory Show memory usage
#' @keywords internal
update_progress <- function(progress, message = "", show_memory = TRUE) {
  if (show_memory) {
    current_mem <- memory_usage()
    mem_diff <- current_mem - progress$memory_start
    message <- sprintf(
      "%s | Memory: %.1f MB (Δ: %+.1f MB)",
      message, current_mem, mem_diff
    )
  }
  
  progress$bar$tick(tokens = list(message = message))
}

#' Get Current Memory Usage
#'
#' @return Memory usage in MB
#' @keywords internal
memory_usage <- function() {
  gc_stats <- gc(full = TRUE, verbose = FALSE)
  sum(gc_stats[, "used"]) * 8 / 1024^2
}

#' Process in Parallel
#'
#' @param X Data to process
#' @param FUN Function to apply
#' @param num_cores Number of cores
#' @param chunk_size Size of chunks
#' @return Results list
#' @importFrom parallel mclapply detectCores
#' @keywords internal
parallel_process <- function(X, FUN, num_cores = 1, chunk_size = NULL) {
  if (num_cores <= 1) {
    return(lapply(X, FUN))
  }
  
  # Set up parallel processing
  num_cores <- min(num_cores, parallel::detectCores())
  
  # Create chunks
  if (is.null(chunk_size)) {
    chunk_size <- ceiling(length(X) / num_cores)
  }
  chunks <- split(X, ceiling(seq_along(X) / chunk_size))
  
  # Process in parallel
  results <- parallel::mclapply(
    chunks,
    function(chunk) lapply(chunk, FUN),
    mc.cores = num_cores
  )
  
  # Combine results
  unlist(results, recursive = FALSE)
}

#' @importFrom methods as
NULL

#' Progress Reporter Constructor
#'
#' @description
#' Creates a progress reporting function with optional timing
#'
#' @param verbose Logical: show messages
#' @param with_timing Logical: include timestamps
#' @return Function for reporting progress
#' @keywords internal
create_reporter <- function(verbose = TRUE, with_timing = TRUE) {
  # Create closure for consistent formatting
  time_format <- if (with_timing) {
    function() sprintf("[%s] ", format(Sys.time(), "%H:%M:%S"))
  } else {
    function() ""
  }
  
  # Return reporter function
  function(message, detail = NULL) {
    if (verbose) {
      cat(sprintf("%s%s\n", time_format(), message))
      if (!is.null(detail)) cat(sprintf("  %s\n", detail))
    }
  }
}

#' Input Validation
#'
#' @description
#' Validates input parameters for CreateCSFNMFobject
#'
#' @param ref_matrix Reference matrix
#' @param ref_celltype Reference cell types
#' @param data_matrix Query matrix
#' @param params List of additional parameters
#' @return TRUE if valid, stops with error otherwise
#' @keywords internal
validate_inputs <- function(ref_matrix, ref_celltype, data_matrix, params) {
  # Check matrix classes
  if (!inherits(ref_matrix, c("matrix", "dgCMatrix"))) {
    stop("ref_matrix must be a matrix or dgCMatrix", call. = FALSE)
  }
  if (!inherits(data_matrix, c("matrix", "dgCMatrix"))) {
    stop("data_matrix must be a matrix or dgCMatrix", call. = FALSE)
  }
  
  # Check dimensions
  if (length(ref_celltype) != ncol(ref_matrix)) {
    stop(sprintf(
      "Length of ref_celltype (%d) must match columns in ref_matrix (%d)",
      length(ref_celltype), ncol(ref_matrix)
    ), call. = FALSE)
  }
  
  # Check parameter values
  if (!params$scale_by %in% c("cells", "genes")) {
    stop('scale_by must be either "cells" or "genes"', call. = FALSE)
  }
  
  # Check cores
  if (params$num_cores > 0) {
    max_cores <- parallel::detectCores()
    if (params$num_cores > max_cores) {
      warning(sprintf(
        "Requested %d cores but only %d available. Using %d cores.",
        params$num_cores, max_cores, max_cores
      ))
      params$num_cores <- max_cores
    }
  }
  
  TRUE
}

#' Safe Matrix Conversion
#'
#' @description
#' Safely converts matrix to sparse format with error handling
#'
#' @param matrix Input matrix
#' @param verbose Logical: show warnings
#' @return Sparse matrix
#' @keywords internal
to_sparse <- function(matrix, verbose = FALSE) {
  tryCatch({
    as(matrix, "CsparseMatrix")
  }, error = function(e) {
    if (verbose) warning("Converting to dense matrix first")
    as(as.matrix(matrix), "CsparseMatrix")
  })
}

#' Memory Usage Reporter
#'
#' @description
#' Reports current memory usage with garbage collection
#'
#' @param threshold Memory threshold in MB to trigger warning
#' @return Named vector of memory statistics
#' @keywords internal
check_memory <- function(threshold = 10000) {
  # Force garbage collection
  gc_stats <- gc(full = TRUE)
  
  # Calculate memory usage
  mem_used <- sum(gc_stats[, "used"]) * 8 / 1024^2  # Convert to MB
  
  # Warning if above threshold
  if (mem_used > threshold) {
    warning(sprintf("High memory usage: %.2f MB", mem_used))
  }
  
  # Return memory stats
  c(
    used = mem_used,
    available = memory.size(max = TRUE) / 1024^2,
    gc_stats = gc_stats
  )
}

#' Cache Manager
#'
#' @description
#' Manages temporary results caching
#'
#' @param max_size Maximum cache size in MB
#' @return Cache environment
#' @keywords internal
create_cache <- function(max_size = 1000) {
  cache <- new.env(parent = emptyenv())
  
  # Add metadata
  cache$size <- 0
  cache$max_size <- max_size
  cache$hits <- 0
  cache$misses <- 0
  
  # Add methods
  cache$set <- function(key, value) {
    if (cache$size >= cache$max_size) {
      warning("Cache full, removing oldest entries")
      rm(list = ls(cache, pattern = "^data_")[1], envir = cache)
    }
    assign(paste0("data_", key), value, envir = cache)
    cache$size <- cache$size + object.size(value) / 1024^2
  }
  
  cache$get <- function(key) {
    key <- paste0("data_", key)
    if (exists(key, envir = cache)) {
      cache$hits <- cache$hits + 1
      get(key, envir = cache)
    } else {
      cache$misses <- cache$misses + 1
      NULL
    }
  }
  
  cache
}