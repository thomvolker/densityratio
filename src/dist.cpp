#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


//' Create a Gram matrix with squared Euclidean distances between
//' observations in the input matrix \code{X} and the input matrix \code{Y}
//' @name distance
//' @param X A numeric input matrix
//' @param Y A numeric input matrix with the same variables as \code{X}
//' @param intercept Logical indicating whether an intercept should be added to
//' the estimation procedure. In this case, the first column is an all-zero
//' column (which will be transformed into an all-ones column in the kernel).
//' @export

// [[Rcpp::export]]
arma::mat distance(const arma::mat& X, const arma::mat& Y, const bool& intercept = false) {
  int nx = X.n_rows;
  int ny = intercept ? Y.n_rows : Y.n_rows - 1;
  int colstart = intercept ? 1 : 0;

  arma::mat XY(nx, ny + 1);
  XY.cols(colstart, ny) -= 2 * X * Y.t();
  XY.cols(colstart, ny).each_col() += sum(X % X, 1);
  XY.cols(colstart, ny).each_row() += sum(Y % Y, 1).t();
  return XY;
}

//' Create gaussian kernel gram matrix from distance matrix
//' @name kernel_gaussian
//' @param dist A numeric distance matrix
//' @param sigma A scalar with the length-scale parameter
//' @export

//[[Rcpp::export]]
arma::mat kernel_gaussian(const arma::mat& dist, const double& sigma) {
  arma::mat KGM = exp(-dist / (2*sigma*sigma));
  return KGM;
}
