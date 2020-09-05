#include "modString.h"
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x*2;
}

//' @title myFunction
//' @description
//' Modify a string in Rcpp.
//' @name myFunction
//' @param x a vector of strings
//' @examples
//' myFunction(x=c('Hello', "C++", 'header', 'files'))
//'
//' @export
// [[Rcpp::export]]
Rcpp::StringVector myFunction(Rcpp::StringVector x) {
  x = modString(x);
  return x;
}

