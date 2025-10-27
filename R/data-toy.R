#' Tiny toy matrices and labels for runnable examples
#'
#' @format
#' - `ref_matrix_toy`: numeric matrix (50 genes x 12 cells), dimnames set.
#' - `qry_matrix_toy`: numeric matrix (50 genes x 8 cells), dimnames set.
#' - `ref_celltype_toy`: character vector of length 12, names match `colnames(ref_matrix_toy)`.
#'
#' @details
#' These are minimal, non-biological toy data with shared gene IDs across
#' reference and query for fast runnable examples and tests.
#'
#' @examples
#' data(ref_matrix_toy, package = "CellMentor")
#' data(qry_matrix_toy, package = "CellMentor")
#' data(ref_celltype_toy, package = "CellMentor")
#' dim(ref_matrix_toy); dim(qry_matrix_toy)
#' head(ref_celltype_toy)
"ref_matrix_toy"

#' @rdname ref_matrix_toy
"qry_matrix_toy"

#' @rdname ref_matrix_toy
"ref_celltype_toy"


#' Tiny prebuilt CSFNMF object for accessors (optional)
#'
#' @format An object of class `csfnmf` built on the toy matrices.
#' @details Only provided to make accessor examples immediate. For real analyses,
#' construct objects from your own data.
#' @examples
#' \donttest{
#' data(obj_toy, package = "CellMentor")
#' inherits(obj_toy, "csfnmf")
#' }
"obj_toy"