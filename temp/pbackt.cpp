#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)
//
// Learn more about how to use Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//
// and browse examples of code using Rcpp at:
//
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector pbackt(NumericVector x) {
  for (int i = 0; i < 10000; i++){
    if (i > 0){
    x.push_back(i);
    } else {
      x[0] = i;
    }
}
  return(x);
}

// [[Rcpp::export]]
NumericVector pbackt2(NumericVector x) {
  for (int i = 0; i < 10000; i++){
    x[i] = i;
  }
  return(x);
}

