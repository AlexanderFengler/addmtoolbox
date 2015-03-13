#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector defaultt(int x = 1, int y = 1,int z = 1) {
  for (int i = 0; i < 1000; i++){
    x++;
    y++;
    z++;
  }

  NumericVector a(3);
  a[0] = x;
  a[1] = y;
  a[2] = z;

  return(a);
}
