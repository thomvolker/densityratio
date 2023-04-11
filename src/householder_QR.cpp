#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double Eucl_Norm(arma::vec X) {
  double squared = 0;
  int k = X.size();
  for (int i = 0; i < k; ++i) {
    double val = X(i);
    squared += val*val;
  }
  return std::sqrt(squared);
}


//' Perform QR decomposition of the input matrix \code{X} through householder
//' transformations, similar to Matlab's QR decomposition.
//'
//' @param X A numeric matrix.
// [[Rcpp::export]]
List householder_QR(arma::mat X) {

  int m = X.n_rows;
  int n = X.n_cols;
  arma::mat Q = arma::mat(m, m, fill::eye);

  int min_dim = std::min(m - 1, n);

  double s, ViTVi;
  arma::vec Vi;
  arma::mat Hi, Qi;

  for (int i = 0; i < min_dim; ++i) {
    arma::mat Xim = X.rows(i, m - 1);
    arma::vec Xi = Xim.col(i);
    if (Xi(0) >= 0) {
      s = 1;
    } else {
      s = -1;
    }
    Vi = zeros(m - i);
    Vi[0] = s * Eucl_Norm(Xi);
    Vi = Vi + Xi;
    ViTVi = as_scalar(Vi.t() * Vi);

    Hi = arma::mat(m - i, m - i, fill::eye) - (Vi * Vi.t()) * (2 / ViTVi);

    if (i == 0) {
      Qi = Hi;
    } else {
      Qi = zeros<mat>(m, m);
      Qi.submat(0, 0, i - 1, i - 1) = arma::mat(i, i, fill::eye);
      Qi.submat(i, i, m - 1, m - 1) = Hi;
    }

    X = Qi * X;
    Q = Q * Qi;

  }
  List out = List::create(Named("Q") = Q, Named("R") = X);
  return out;
}
