#include <Rcpp.h>
using namespace Rcpp;

NumericVector utilities_timesTwo(NumericVector x) {
  return x * 4;
}

