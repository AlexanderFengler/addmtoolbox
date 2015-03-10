#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix fixmodt(NumericVector x, int size, Function fixation_model) {
  NumericMatrix y = fixation_model(x,size);
  return y;
}
