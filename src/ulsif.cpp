#include <RcppArmadillo.h>
#include <omp.h>
using namespace Rcpp;
using namespace arma;

arma::vec compute_alpha(arma::mat Hhat, arma::vec hhat, arma::mat diagmat, double lambda) {
  arma::vec alpha = solve(Hhat + lambda * diagmat, hhat);
  return alpha;
}

int set_threads(int nthreads) {
  int max_threads = omp_get_max_threads();
  if (nthreads == 0) {
    nthreads = max_threads;
  }
  if (nthreads > max_threads || nthreads == -1) {
    nthreads = max_threads;
    std::string warn = "'nthreads' must specify the number of threads to use for parallel processing; it is set to the maximum number of threads available (" + std::to_string(nthreads) + ")";
    Rf_warningcall(R_NilValue, warn.c_str());
  }
  return nthreads;
}

//[[Rcpp::export]]
arma::mat compute_ulsif(arma::mat &Hhat, arma::vec &hhat, arma::vec &lambda, bool parallel, int nthreads) {
  int nlambda = lambda.size();
  int p = Hhat.n_cols;
  arma::mat diagmat = arma::eye<arma::mat>(p, p);
  arma::mat alpha = arma::zeros<arma::mat>(p, nlambda);

  if (parallel) {
    nthreads = set_threads(nthreads);
  }
  #pragma omp parallel for num_threads(nthreads) if (parallel)
  for(int i = 0; i < nlambda; i++) {
    alpha.col(i) = compute_alpha(Hhat, hhat, diagmat, lambda(i));
  }
  return alpha;
}

