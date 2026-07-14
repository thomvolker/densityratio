#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat distance(const arma::mat& X, const arma::mat& Y, const bool& intercept = false) {
  int nx = X.n_rows;
  int ny = intercept ? Y.n_rows : Y.n_rows - 1;
  int colstart = intercept ? 1 : 0;

  arma::mat XY(nx, ny + 1, arma::fill::zeros);
  XY.cols(colstart, ny) -= 2 * X * Y.t();
  XY.cols(colstart, ny).each_col() += sum(X % X, 1);
  XY.cols(colstart, ny).each_row() += sum(Y % Y, 1).t();
  return XY;
}

//[[Rcpp::export]]
arma::mat kernel_gaussian(const arma::mat& dist, double sigma) {
  arma::mat KGM = exp(-dist / (2*sigma*sigma));
  return KGM;
}
