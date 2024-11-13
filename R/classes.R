#' @importFrom methods setClass new
#' @importFrom Matrix sparseMatrix
NULL

#' RefDataList Class
#' 
#' @description
#' Class for storing paired reference and query expression matrices
#' 
#' @slot ref Reference matrix (sparse Matrix) containing gene expression from labeled cells
#' @slot data Query matrix (sparse Matrix) containing gene expression from unlabeled cells
#' @export
setClass("RefDataList", 
         slots = c(
           ref = "Matrix", 
           data = "Matrix"
         ))

#' Annotation Class
#' 
#' @description
#' Class for storing cell type annotations and their numeric encodings
#' 
#' @slot ref_celltype Vector of cell type labels
#' @slot ref_celltype.code Numeric vector of encoded cell types
#' @export
setClass("Annotation", 
         slots = c(
           ref_celltype = "vector", 
           ref_celltype.code = "vector"
         ))

#' Rank Class
#' 
#' @description
#' Class for storing rank information
#' 
#' @slot k Numeric rank value
#' @export
setClass("rank", 
         slots = c(
           k = "numeric"
         ))

#' Helper Matrices Class
#' 
#' @description
#' Class for storing intermediate matrices used in NMF computation
#' 
#' @slot A First matrix
#' @slot B Second matrix
#' @slot P P matrix
#' @slot M M matrix
#' @slot N N matrix
#' @slot BP BP matrix
#' @slot BP_posneg List of positive/negative BP values
#' @slot Hconst List of H constants
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
#' Class for storing training results and parameters
#' 
#' @slot count.matrices RefDataList object
#' @slot annotation Data frame of annotations
#' @slot test.annotation Data frame of test annotations
#' @slot rank Rank object
#' @slot H H matrix
#' @slot W W matrix
#' @slot constants Helper matrices
#' @slot hyper_para List of hyperparameters
#' @slot init_method Character indicating initialization method
#' @slot max.iter Numeric maximum iterations
#' @slot loss Numeric loss value
#' @slot accuracy Numeric accuracy value
#' @export
setClass("traincsfnmf",
         slots = c(
           count.matrices = "RefDataList",
           annotation = "data.frame",
           test.annotation = "data.frame",
           rank = "rank",
           H = "Matrix",
           W = "Matrix",
           constants = "helpmat",
           hyper_para = "list",
           init_method = "character",
           max.iter = "numeric",
           loss = "numeric",
           accuracy = "numeric"
         ))

#' CSFNMF Class
#' 
#' @description
#' Main class for Constrained Supervised Factorization NMF analysis
#' 
#' @slot count.matrices Current working matrices
#' @slot origin.matrices Original unprocessed matrices
#' @slot norm.matrices Normalized matrices
#' @slot scale.matrices Scaled matrices
#' @slot f.select.matrices Feature-selected matrices
#' @slot annotation Data frame of annotations
#' @slot train_object Trained model parameters
#' @slot H Factor loading matrix for cells
#' @slot W Factor loading matrix for genes
#' @export
setClass("csfnmf",
         slots = c(
           count.matrices = "RefDataList",
           origin.matrices = "RefDataList",
           norm.matrices = "RefDataList",
           scale.matrices = "RefDataList",
           f.select.matrices = "RefDataList",
           annotation = "data.frame",
           train_object = "traincsfnmf",
           H = "Matrix",
           W = "Matrix"
         ))

#' Extended Rank Class
#' @description Class for storing rank information with additional parameters
#' @export
setClass("ExtendedRank",
         contains = "rank",
         slots = c(k = "numeric",
                   alpha = "numeric",
                   beta = "numeric",
                   scale_vectors = "list"))

#' Extended Training CSFNMF Class
#' 
#' @description Class for storing training results with updates
#' @exportClass update_traincsfnmf
setClass("update_traincsfnmf", 
         contains = "traincsfnmf",
         slots = c(updates = "list"))