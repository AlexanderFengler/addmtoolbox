#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
float vect(NumericVector x) {
  float x1 = x[0];
  float x2 = x[1];
  float x3 = x[2];
  for (int i = 0; i < 10000; i++){
    x1++;
    x2++;
    x3++;
  }
  return x1;
}


// [[Rcpp::export]]
float vect2(float x1 = 0, float  x2 = 0, float x3 = 0) {
  for (int i = 0; i < 10000; i++){
    x1++;
    x2++;
    x3++;
  }
  return x1;
}
