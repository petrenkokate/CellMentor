#' @importFrom methods setClass new
#' @importFrom Matrix Matrix
NULL

#' RefDataList Class
#' 
#' @description
#' Class for storing paired reference and query expression matrices
#' 
#' @slot ref Reference matrix
#' @slot data Query matrix
#' @export
setClass("RefDataList", 
         slots = c(
           ref = "Matrix", 
           data = "Matrix"
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
#' @slot matrices Current working matrices (RefDataList)
#' @slot annotation Training data annotations
#' @slot test_annotation Test data annotations
#' @slot rank Numeric rank value
#' @slot H H matrix
#' @slot W W matrix
#' @slot constants Helper matrices
#' @slot parameters List of parameters (hyperparameters, initialization method, etc.)
#' @slot results List of results (loss, accuracy)
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
#' Main class for Constrained Supervised Factorization NMF analysis
#' 
#' @slot matrices Current working matrices (RefDataList)
#' @slot annotation Data frame of annotations
#' @slot train_object Training results and parameters
#' @slot H Final H matrix
#' @slot W Final W matrix
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
#' Class for storing training results with updates
#' 
#' @slot base Base training results
#' @slot updates List of update results
#' @export
setClass("update_traincsfnmf", 
         contains = "traincsfnmf",
         slots = c(updates = "list"))