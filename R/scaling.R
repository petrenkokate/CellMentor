#' @importFrom Matrix rowSums colSums
NULL

#' Scale Matrix by RMSD
#'
#' @description
#' Scales matrix by root mean square deviation
#'
#' @param matrix Input matrix
#' @param dimension "row" or "column"
#' @param eps Small value to prevent division by zero
#' @return Scaled matrix
#' @keywords internal
#' @noRd
scale_rmsd <- function(matrix, dimension = c("row", "column"), eps = 1e-10) {
  dimension <- match.arg(dimension)

  # Determine scaling parameters
  if (dimension == "column") {
    sums <- Matrix::colSums(matrix * matrix)
    n <- nrow(matrix) - 1L
    margin <- 2L
  } else {
    sums <- Matrix::rowSums(matrix * matrix)
    n <- ncol(matrix) - 1L
    margin <- 1L
  }

  # Calculate scaling factors
  root_mean_sq <- sqrt(sums / n)
  root_mean_sq[root_mean_sq < eps] <- eps

  # Scale matrix
  scaled <- sweep(matrix, margin, root_mean_sq, "/", check.margin = FALSE)

  # Ensure sparse format
  as(scaled, "CsparseMatrix")
}

#' Scale Matrix by Cells
#'
#' @description
#' Scales matrix by cells using RMSD
#'
#' @param matrix Input matrix
#' @param eps Small value to prevent division by zero
#' @return Scaled matrix
#' @keywords internal
#' @noRd
scale_by_cells <- function(matrix, eps = 1e-10) {
  scale_rmsd(matrix, "column", eps)
}

#' Scale Matrix by Genes
#'
#' @description
#' Scales matrix by genes using RMSD
#'
#' @param matrix Input matrix
#' @param eps Small value to prevent division by zero
#' @return Scaled matrix
#' @keywords internal
#' @noRd
scale_by_genes <- function(matrix, eps = 1e-10) {
  scale_rmsd(matrix, "row", eps)
}

#' Normalize Matrix
#'
#' @description
#' Performs log normalization on expression matrix
#'
#' @param matrix Input matrix
#' @param scale.factor Scale factor for normalization
#' @return Normalized matrix
#' @keywords internal
#' @noRd
normalize_matrix <- function(matrix, scale.factor = 1e4) {
  # Calculate size factors
  size_factors <- Matrix::colSums(matrix)
  size_factors[size_factors == 0] <- 1

  # Normalize and log transform efficiently
  normalized <- 
    sweep(matrix, 2, size_factors, "/", check.margin = FALSE) * scale.factor
  log1p(normalized)
}
