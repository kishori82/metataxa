#' @title hello
#' @description
#' Says hello.
#' @name hello
#' @param txt something to say hello to
#'
#' @export

hello <- function(txt = "world") {
  cat("Hello, ", txt, "\n")
}


#' The first test function
#'
#' @param a A factor
#' @param b Another factor
#'
#' @return factor
#' @export
#'


fbind <- function(a, b) {
  factor(c(as.character(a), as.character(b)))
}


#' @title delta_star
#' @description
#' @name delta_star
#' @param file_name file name with taxa
#'
#' @export

delta_star <- function(file_name = NULL) {
  compute_tax_distance(file_name)
  dir()
  return(read.csv(file_name))
}

#' @useDynLib metataxa
#' @importFrom Rcpp sourceCpp
NULL

