#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int sizetest(NumericVector x) {
  int y = x.size();
  return y;
}
