#' @title vcovConley-package
#' 
#' @description Functions to calculate conley standard errors.
#' 
#' @author Nicholas Potter
#' @docType package
#' @name vcovConley
#' @aliases vcovConley
#' @useDynLib vcovConley, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' FUnction to unload compiled code
.onUnload <- function (libpath) {
  library.dynam.unload("vcovConley", libpath)
}