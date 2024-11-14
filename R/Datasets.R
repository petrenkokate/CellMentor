#' @importFrom methods setClass new
#' @importFrom Matrix Matrix
NULL

#' RefDataList Class
#' 
#' @description
#' Class for storing paired reference and query expression matrices. The reference matrix contains
#' gene expression data from labeled cells, while the data matrix contains expression data from
#' unlabeled cells that need to be classified.
#' 
#' @slot ref Reference matrix (genes × cells) containing expression data from labeled cells
#' @slot data Query matrix (genes × cells) containing expression data from unlabeled cells
#' @export
setClass("RefDataList", 
         slots = c(
           ref = "Matrix", 
           data = "Matrix"
         ))

#' Helper Matrices Class
#' 
#' @description
#' Class for storing intermediate matrices used in NMF computation. These matrices
#' are used to optimize the factorization process and enforce constraints.
#' 
#' @slot A Matrix for within-class scatter computation
#' @slot B Matrix for between-class scatter computation
#' @slot P Projection matrix for class relationships
#' @slot M Regularization matrix for W updates
#' @slot N Normalization matrix for class weights
#' @slot BP Combined B and P matrices for efficiency
#' @slot BP_posneg List containing positive and negative components of BP matrix
#' @slot Hconst List of constants used in H matrix updates
#' @export
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
         ))

#' TrainCSFNMF Class
#' 
#' @description
#' Class for storing training results and parameters. This class manages the training
#' process where reference data is split into training and validation sets to optimize
#' the NMF model.
#' 
#' @slot matrices RefDataList containing:
#'        ref: Training portion of reference data
#'        data: Validation portion of reference data
#' @slot annotation Data frame containing cell type labels for training data
#' @slot test_annotation Data frame containing cell type labels for validation data
#' @slot rank Numeric value indicating the factorization rank (number of factors)
#' @slot H Factor matrix (rank × cells) representing cell type signatures
#' @slot W Factor matrix (genes × rank) representing gene weights
#' @slot constants Helper matrices for computation optimization
#' @slot parameters List of model parameters including:
#'        - rank: factorization rank
#'        - init_method: initialization method
#'        - alpha: within-class scatter weight
#'        - beta: between-class scatter weight
#'        - gamma: sparsity parameter
#'        - delta: orthogonality parameter
#' @slot results List containing training outcomes:
#'        - loss: optimization loss values
#'        - accuracy: classification accuracy
#' @export
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
         ))

#' CSFNMF Class
#' 
#' @description
#' Main class for Constrained Supervised Factorization NMF analysis. This class
#' orchestrates the entire analysis workflow, from initial data organization to
#' final cell type prediction.
#' 
#' @slot matrices RefDataList containing:
#'        ref: Complete reference dataset with known cell types
#'        data: Query dataset to be classified
#' @slot annotation Data frame containing cell type labels for the complete reference dataset
#' @slot train_object TrainCSFNMF object containing training results and model parameters
#' @slot H Final factor matrix for all cells (combined reference and query)
#' @slot W Final gene weight matrix
#' @export
setClass("csfnmf",
         slots = c(
           matrices = "RefDataList",
           annotation = "data.frame",
           train_object = "traincsfnmf",
           H = "Matrix",
           W = "Matrix"
         ))

#' Training Results Class
#' 
#' @description
#' Class for storing training results with updates. This class extends TrainCSFNMF
#' to track model performance across multiple training iterations and parameter updates.
#' 
#' @slot base Inherited slots from TrainCSFNMF
#' @slot updates List containing training history:
#'        - WH_list: Factor matrices for each update
#'        - loss_list: Loss values across updates
#'        - accuracy_list: Accuracy measures for each update
#' @export
setClass("update_traincsfnmf", 
         contains = "traincsfnmf",
         slots = c(updates = "list"))