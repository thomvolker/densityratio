#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Create a Gram matrix with squared Euclidean distances between all
//' observations in the input matrix
//'
//' @param x A numeric input matrix

// [[Rcpp::export]]
arma::mat distX(arma::mat X) {
  double dev, dist;
  int nr = X.n_rows;
  int nc = X.n_cols;
  arma::mat out = arma::zeros<arma::mat>(nr, nr);

  for(int i = 0; i < nr; i++) {
    for(int j = i+1; j < nr; j++) {
      dist = 0;
      for(int c = 0; c < nc; c++) {
        dev = X(i,c) - X(j,c);
        dist += dev*dev;
      }
      out(i,j) = dist;
    }
  }
  return out + out.t();
}
