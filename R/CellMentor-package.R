#' CellMentor: Cell-Type Aware Dimensionality Reduction for Single-Cell RNA-Seq
#'
#' CellMentor is a supervised dimensionality reduction method based on
#' non-negative matrix factorization (NMF) that integrates cell type labels
#' directly into its optimization objective. By minimizing variation within
#' known populations while maximizing distinctions between types, CellMentor
#' produces low-dimensional embeddings optimized for cell type identification
#' in single-cell RNA sequencing analysis.
#'
#' @section Key Features:
#' \itemize{
#'   \item \strong{Supervised NMF Framework:} Incorporates labels via discriminative constraints
#'   \item \strong{Superior Cell Type Separation:} Maximally separable embeddings
#'   \item \strong{Robust Batch Handling:} Preserves biology while mitigating technical effects
#'   \item \strong{Rare Population Detection:} Sensitive to low-frequency types
#'   \item \strong{Automated Parameter Optimization:} Built-in hyperparameter tuning
#' }
#'
#' @section Two-Phase Workflow:
#' \enumerate{
#'   \item \strong{Decomposition (Training):} Learn W (genes × K) and H (K × cells)
#'   \item \strong{Projection (Inference):} Project queries with non-negative least squares
#' }
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{CellMentor}} — supervised NMF / hyperparameter search
#'   \item \code{\link{project_data}} — project queries using learned W
#'   \item \code{\link{CreateCSFNMFobject}} — initialize a CellMentor object
#' }
#'
#' @section Getting Help:
#' \itemize{
#'   \item Package docs: \code{help(package = "CellMentor")}
#'   \item GitHub: \url{https://github.com/petrenkokate/CellMentor}
#'   \item Issues: \url{https://github.com/petrenkokate/CellMentor/issues}
#' }
#'
#' @author
#' Or Hevdeli (equal contribution) \cr
#' Ekaterina Petrenko (equal contribution) \email{petrenko.kate@icloud.com} \cr
#' Dvir Aran (corresponding author) \email{dviraran@technion.ac.il}
#'
#' @references
#' Hevdeli, O., Petrenko, E., & Aran, D. (2025). CellMentor: Cell-Type Aware
#' Dimensionality Reduction for Single-cell RNA-Sequencing Data. \emph{bioRxiv}.
#' \doi{10.1101/2025.06.17.660094}
#'
#' @name CellMentor
#' @docType package
#' @aliases CellMentor-package
"_PACKAGE"
