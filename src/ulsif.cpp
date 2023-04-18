#include <RcppArmadillo.h>
#include <omp.h>
using namespace Rcpp;
using namespace arma;

arma::vec compute_alpha(arma::mat Hhat, arma::vec hhat, arma::mat diagmat, double lambda) {
  arma::vec alpha = solve(Hhat + lambda * diagmat, hhat);
  return alpha;
}

//[[Rcpp::export]]
arma::mat Cpp_ulsif(arma::mat &Hhat, arma::vec &hhat, arma::vec &lambda, bool parallel, int nthreads) {
  int nlambda = lambda.size();
  int p = Hhat.n_cols;
  arma::mat diagmat = arma::eye<arma::mat>(p, p);
  arma::mat alpha = arma::zeros<arma::mat>(p, nlambda);

  if (parallel) {
    int max_threads = omp_get_max_threads();
    if (nthreads == 0 || nthreads > max_threads) {
      nthreads = omp_get_max_threads();
    }
  }
  #pragma omp parallel for num_threads(nthreads) if (parallel)
  for(int i = 0; i < nlambda; i++) {
    alpha.col(i) = compute_alpha(Hhat, hhat, diagmat, lambda(i));
  }
  return alpha;
}

