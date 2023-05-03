#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


//' Create a Gram matrix with squared Euclidean distances between
//' observations in the input matrix \code{X} and the input matrix \code{Y}
//' @param X A numeric input matrix
//' @param Y A numeric input matrix with the same variables as \code{X}
//' @param symmetric A logical indicating whether X and Y are the same

// [[Rcpp::export]]
arma::mat distance(arma::mat X, arma::mat Y, bool symmetric = false) {
  double dev, dist;
  int nrx = X.n_rows;
  int nry = Y.n_rows;
  int nc = X.n_cols;

  arma::mat out(nrx, nry);

  if (symmetric) {
    for (int i = 0; i < nrx; i++) {
      for(int j = i + 1; j < nry; j++) {
        dist = 0;
        for(int c = 0; c < nc; c++) {
          dev = X(i,c) - Y(j,c);
          dist += dev*dev;
        }
        out(i,j) = dist;
      }
    }
    out = out + out.t();
  } else {
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
  }
  return out;
}

arma::vec make_sigma(arma::mat dist, arma::vec sigma) {

  arma::vec distvec = arma::vectorise(dist);
  arma::vec P, Q, s;

  if (sigma(0) == 0) {
    s = arma::median(distvec);
  } else if (sigma(0) == -1) {
    P = {0.25, 0.3, 0.35, 0.4, 0.45, 0.475, 0.5, 0.525, 0.55, 0.6, 0.65, 0.7, 0.75};
    Q = quantile(distvec, P);
    s = unique(Q);
  } else {
    s = sigma;
  }
  return s;
}

arma::mat kernel_gaussian(arma::mat dist, double sigma) {
  arma::mat KGM = exp(-dist / (2*sigma*sigma));
  return KGM;
}




