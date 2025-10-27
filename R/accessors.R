#' Access the W matrix
#'
#' Retrieve the W matrix (gene loadings) from a CellMentor object.
#'
#' @param x A CellMentor object (e.g., `traincsfnmf` or `csfnmf`).
#' @return A numeric matrix representing the W (gene loadings) matrix.
#' @examples
#' data(obj_toy, package = "CellMentor")
#' W(obj_toy)
#' @export
W <- function(x) x@W


#' Access the H matrix
#'
#' Retrieve the H matrix (cell embeddings) from a CellMentor object.
#'
#' @param x A CellMentor object.
#' @return A numeric matrix representing the H (cell embeddings) matrix.
#' @examples
#' data(obj_toy, package = "CellMentor")
#' H(obj_toy)
#' @export
H <- function(x) x@H


#' Access the rank of a CellMentor object
#'
#' Retrieve the factorization rank used during CSFNMF training.
#'
#' @param x A CellMentor object.
#' @return A numeric value representing the selected rank.
#' @examples
#' \dontrun{
#' # Access the rank of the model
#' cm_rank(cs_obj)
#' }
#' @export
cm_rank <- function(x) x@parameters$rank


#' Access cell annotations from a CellMentor object
#'
#' Retrieve cell type or metadata annotations.
#'
#' @param x A CellMentor object.
#' @return A data frame containing cell annotations.
#' @examples
#' data(obj_toy, package = "CellMentor")
#' cm_annotation(obj_toy)
#' @export
cm_annotation <- function(x) x@annotation


#' Access the RefDataList from a CSFNMF object
#'
#' Retrieve the RefDataList structure containing both reference and 
#' query matrices.
#'
#' @param x A `csfnmf` or `traincsfnmf` object.
#' @return A RefDataList object containing the reference and query matrices.

#' @examples
#' data(obj_toy, package = "CellMentor")
#' matrices(obj_toy)
#' @export
matrices <- function(x) x@matrices


#' Access the query (data) matrix from a RefDataList object
#'
#' Retrieve the single-cell expression matrix for the query dataset.
#'
#' @param x A RefDataList object.
#' @return A sparse `Matrix::Matrix` object representing query data.
#' @examples
#' data(obj_toy, package = "CellMentor")
#' data_matrix(matrices(obj_toy))
#' @export
data_matrix <- function(x) x@data


#' Access the reference matrix from a RefDataList object
#'
#' Retrieve the single-cell expression matrix for the reference dataset.
#'
#' @param x A RefDataList object.
#' @return A sparse `Matrix::Matrix` object representing reference data.
#' @examples
#' data(obj_toy, package = "CellMentor")
#' ref_matrix(matrices(obj_toy))
#' @export
ref_matrix <- function(x) x@ref
