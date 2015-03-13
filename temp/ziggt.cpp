// [[Rcpp::depends(RcppZiggurat)]]
#include <Rcpp.h>
#include <algorithm>
#include <stdlib.h>
#include <Ziggurat.h>
using namespace Rcpp;
static Ziggurat::Ziggurat::Ziggurat zigg;

// [[Rcpp::export]]
NumericVector ziggt() {
  NumericVector x  = zigg.gsl(5);
  return x;
}
