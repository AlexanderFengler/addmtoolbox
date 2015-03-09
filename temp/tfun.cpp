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
NumericVector tfun(int nr_attributes,
                   int nr_items = 2,
                   NumericVector update = 3,
                   float theta = 0.5,
                   float gamma = 0.5,
                   float drift = 0.002,
                   int x = 0) {

  NumericMatrix updatemat(nr_attributes,2);
  NumericVector drifts(2*nr_attributes);

  for(int i = 0; i < nr_items; i++){
    for (int j = 0; j < nr_attributes; j++){
      updatemat(j,i) = update[((nr_attributes*i)+j)];
    }
  }

  for (int i = 0; i < nr_items; i++){
    for (int j = 0; j < nr_attributes; j++){
      for (int k = 0; k <nr_attributes; k++){
        if (j == k){
          if (i == 0){
            drifts[((nr_attributes*i)+j)] += updatemat(j,i) - theta*updatemat(j,1);
          } else {
            drifts[((nr_attributes*i)+j)] += theta*updatemat(j,0) - updatemat(j,i);
          }
        } else {
          if (i == 0){
            drifts[((nr_attributes*i)+j)] += gamma*(updatemat(j,i) - theta*updatemat(j,1));
          } else {
            drifts[((nr_attributes*i)+j)] += gamma*(theta*updatemat(j,0) - updatemat(j,i));
          }
        }
      }
    }
  }

  for (int i = 0; i < drifts.size(); i++){
    drifts[i] = drifts[i]*drift;
  }

 return(drifts);
}
