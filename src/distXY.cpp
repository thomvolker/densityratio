#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


//' Create a Gram matrix with squared Euclidean distances between
//' observations in the input matrix \code{X} and the input matrix \code{Y}
//' @param X A numeric input matrix
//' @param Y A numeric input matrix with the same variables as x

// [[Rcpp::export]]
arma::mat distXY(arma::mat X, arma::mat Y) {
  double dev, dist;
  int nrx = X.n_rows;
  int nry = Y.n_rows;
  int nc = X.n_cols;

  arma::mat out(nrx, nry);

  for(int i = 0; i < nrx; i++) {
    for(int j = 0; j < nry; j++) {
      dist = 0;
      for (int c = 0; c < nc; c++) {
        dev = X(i,c) - Y(j,c);
        dist += dev*dev;
      }
      out(i,j) = dist;
    }
  }
  return out;
}




