#' @importFrom methods setClass new
#' @importFrom Matrix Matrix
NULL

#' RefDataList Class
#' 
#' @description
#' An S4 class for storing paired reference and query expression matrices
#' used in CSFNMF analysis. This class ensures consistent handling of
#' sparse matrix data throughout the factorization pipeline.
#' 
#' @details
#' The RefDataList class is designed to maintain both reference (training)
#' and query (test) single-cell RNA-seq expression matrices in sparse format.
#' Matrices are stored as Matrix objects to efficiently handle the sparsity
#' typical of scRNA-seq data.
#' 
#' Structure:
#' \itemize{
#'   \item Both matrices must have the same gene set (rows)
#'   \item Genes are rows, cells are columns
#'   \item Data is typically log-normalized expression values
#'   \item Sparse Matrix format (dgCMatrix) for memory efficiency
#' }
#' 
#' @slot ref Matrix object containing the reference expression matrix.
#'   Dimensions: genes (rows) × reference cells (columns).
#'   Contains normalized expression values from the reference dataset
#'   with known cell type annotations.
#' @slot data Matrix object containing the query/data expression matrix.
#'   Dimensions: genes (rows) × query cells (columns).
#'   Contains normalized expression values from the dataset to be
#'   annotated or analyzed.
#' @noRd
setClass("RefDataList",
  slots = c(
    ref = "Matrix",
    data = "Matrix"
  )
)

#' Helper Matrices Class
#' 
#' @description
#' An S4 class for storing intermediate matrices and constants used during
#' CSFNMF optimization. These matrices are pre-computed to improve
#' computational efficiency during the iterative update steps.
#' 
#' @details
#' The helpmat class contains auxiliary matrices that remain constant
#' throughout the optimization iterations. Pre-computing these values
#' significantly reduces computational overhead in the multiplicative
#' update rules.
#' 
#' Mathematical Context:
#' \itemize{
#'   \item A, B: Constraint matrices for supervised factorization
#'   \item P, M, N: Scatter matrices for within/between class constraints
#'   \item BP: Product matrix used in H updates
#'   \item BP_posneg: Positive and negative components of BP for multiplicative updates
#'   \item Hconst: Pre-computed constants for H matrix updates
#' }
#' 
#' @slot A Matrix object. First constraint matrix used in the factorization.
#' @slot B Matrix object. Second constraint matrix used in the factorization.
#' @slot P Matrix object. Pairwise constraint matrix.
#' @slot M Matrix object. Within-class scatter matrix for supervised constraints.
#' @slot N Matrix object. Between-class scatter matrix for supervised constraints.
#' @slot BP Matrix object. Pre-computed product of B and P matrices.
#' @slot BP_posneg List containing positive and negative components of BP.
#'   Structure: list(pos = positive_part, neg = negative_part).
#'   Used in multiplicative update rules to ensure non-negativity.
#' @slot Hconst List of constants for H matrix updates.
#'   Pre-computed values that remain constant during optimization iterations.
#' @noRd
setClass("helpmat",
  slots = c(
    A = "Matrix",
    B = "Matrix",
    P = "Matrix",
    M = "Matrix",
    N = "Matrix",
    BP = "Matrix",
    BP_posneg = "list",
    Hconst = "list"
  )
)

#' TrainCSFNMF Class
#' 
#' @description
#' An S4 class for storing training results, parameters, and intermediate
#' matrices during CSFNMF model fitting. This class encapsulates all
#' information needed for the optimization process and final model evaluation.
#' 
#' @details
#' The traincsfnmf class represents a complete CSFNMF model after training.
#' It contains the factorization results (W and H matrices), training
#' parameters, optimization results, and the data used for training.
#' 
#' Workflow Integration:
#' \itemize{
#'   \item Created by divide_reference_data() from a csfnmf object
#'   \item Updated iteratively by RunCSFNMF() during optimization
#'   \item Contains convergence metrics and final loss values
#'   \item Used for projecting new data via project_data()
#' }
#' 
#' @slot matrices RefDataList object containing the working matrices.
#'   Contains both reference and data matrices used in training.
#'   See ?RefDataList for structure details.
#' @slot annotation data.frame of cell type annotations for training data.
#'   Contains at minimum a 'celltype' column with cell type labels.
#'   Row names correspond to cell identifiers in the reference matrix.
#' @slot test_annotation data.frame of cell type annotations for test data.
#'   Used for evaluating model performance during training.
#'   Same structure as annotation slot.
#' @slot rank Numeric value specifying the factorization rank (k).
#'   Determines the number of latent factors/components in the decomposition.
#'   Typically ranges from 20-100 depending on dataset complexity.
#' @slot H Matrix object containing the learned H matrix (cell embeddings).
#'   Dimensions: k (rank) × cells.
#'   Represents cell embeddings in the low-dimensional latent space.
#' @slot W Matrix object containing the learned W matrix (gene loadings).
#'   Dimensions: genes × k (rank).
#'   Represents gene weights for each latent factor.
#' @slot constants helpmat object containing pre-computed helper matrices.
#'   See ?helpmat for details on the auxiliary matrices used in optimization.
#' @slot parameters List of hyperparameters used in training.
#'   Contains: rank, init_method, alpha, beta, gamma, delta, max_iter.
#'   Controls the optimization behavior and constraints.
#' @slot results List of training results and convergence metrics.
#'   Contains: loss (vector of loss values per iteration),
#'   nmi (Normalized Mutual Information score),
#'   predictions (predicted cell types for evaluation).
#' @noRd
setClass("traincsfnmf",
  slots = c(
    matrices = "RefDataList",
    annotation = "data.frame",
    test_annotation = "data.frame",
    rank = "numeric",
    H = "Matrix",
    W = "Matrix",
    constants = "helpmat",
    parameters = "list",
    results = "list"
  )
)

#' CSFNMF Class
#'
#' @description
#' Main class for Constrained Supervised Factorization NMF analysis
#'
#' @slot matrices Current working matrices (RefDataList)
#' @slot annotation Data frame of annotations
#' @slot train_object Training results and parameters
#' @slot H Final H matrix
#' @slot W Final W matrix
#' @noRd
setClass("csfnmf",
  slots = c(
    matrices = "RefDataList",
    annotation = "data.frame",
    train_object = "traincsfnmf",
    H = "Matrix",
    W = "Matrix"
  )
)

#' Training Results with Updates Class
#' 
#' @description
#' An S4 class extending traincsfnmf to track optimization progress across
#' multiple iterations or updates. This class stores both the current state
#' and historical snapshots of the training process.
#' 
#' @details
#' The update_traincsfnmf class inherits all slots from traincsfnmf and adds
#' an updates slot to maintain a history of model states. This is useful for:
#' \itemize{
#'   \item Monitoring convergence behavior over iterations
#'   \item Comparing different initialization strategies
#'   \item Debugging optimization issues
#'   \item Analyzing parameter sensitivity
#' }
#' 
#' Usage Context:
#' This class is primarily used internally during iterative optimization
#' to maintain checkpoints of the model state. Most users will work with
#' the standard traincsfnmf class returned by RunCSFNMF().
#' 
#' @slot base Inherited from traincsfnmf. See ?traincsfnmf for all inherited slots:
#'   matrices, annotation, test_annotation, rank, H, W, constants,
#'   parameters, results.
#' @slot updates List of historical model states.
#'   Each element is a snapshot containing: iteration number, loss value,
#'   parameter values, and optionally W/H matrices at that iteration.
#'   Used for convergence diagnostics and model comparison.
#'
#' @noRd
setClass("update_traincsfnmf",
  contains = "traincsfnmf",
  slots = c(updates = "list")
)
