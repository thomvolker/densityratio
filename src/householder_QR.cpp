#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


//' Perform QR decomposition of the input matrix \code{X} through householder
//' transformations, similar to Matlab's QR decomposition.
//'
//' @param X A numeric matrix.

// [[Rcpp::export]]
List householder_QR(arma::mat X) {

  int m = X.n_rows;
  int n = X.n_cols;
  arma::mat Q = arma::eye<arma::mat>(m, m);

  int min_dim = std::min(m - 1, n);

  double s, ViTVi;
  arma::vec Xi, Vi;
  arma::mat Hi, Qi;

  for (int i = 0; i < min_dim; ++i) {
    Xi = X.submat(i, i, m - 1, i);
    double s = (Xi(0) >= 0) ? 1.0 : -1.0;

    Vi = arma::zeros<arma::vec>(m - i);
    Vi(0) = s * arma::norm(Xi, 2);
    Vi += Xi;
    ViTVi = arma::dot(Vi, Vi);

    Hi = arma::eye<arma::mat>(m - i, m - i) - (Vi * Vi.t()) * (2 / ViTVi);

    if (i == 0) {
      Qi = Hi;
    } else {
      Qi = zeros<mat>(m, m);
      Qi.submat(0, 0, i - 1, i - 1) = arma::eye<arma::mat>(i, i);
      Qi.submat(i, i, m - 1, m - 1) = Hi;
    }

    X = Qi * X;
    Q *= Qi;

  }
  List out = List::create(Named("Q") = Q, Named("R") = X);
  return out;
}
