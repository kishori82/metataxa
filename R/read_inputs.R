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
#' @useDynLib metataxa
#' @importFrom Rcpp sourceCpp
NULL


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

