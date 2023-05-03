#include <RcppArmadillo.h>
#include <omp.h>
#include "densityratio.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::vec compute_alpha(arma::mat Hhat, const arma::vec& hhat, const double lambda) {
  Hhat.diag() += lambda;
  arma::vec alpha = arma::solve(Hhat, hhat);
  return alpha;
}

//[[Rcpp::export]]
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
arma::mat make_Hhat(arma::mat dist_de, double sigma) {
  int n_de = dist_de.n_rows;
  arma::mat phi_de = kernel_gaussian(dist_de, sigma);
  arma::mat Hhat = phi_de.t() * phi_de / n_de;
  return Hhat;
}

//[[Rcpp::export]]
arma::mat make_hhat(arma::mat dist_nu, double sigma) {
  arma::mat phi_nu = kernel_gaussian(dist_nu, sigma);
  arma::mat hhat = arma::mean(phi_nu, 0);
  return hhat.t();
}

//[[Rcpp::export]]
List compute_ulsif(arma::mat dist_nu, arma::mat dist_de, arma::vec sigma, arma::vec lambda, bool parallel, int nthreads) {

  arma::vec s = make_sigma(dist_nu, sigma);
  double si;
  int nsigma  = s.size();
  int nlambda = lambda.size();
  int nc = dist_nu.n_cols;
  int sig, l;

  arma::mat diagmat = arma::eye<arma::mat>(nc, nc);
  arma::mat Hhat = arma::zeros<arma::mat>(nc, nc);
  arma::mat hhat = arma::zeros<arma::mat>(nc);
  arma::cube alpha(nc, nlambda, nsigma, arma::fill::zeros);

  if (parallel) {
    nthreads = set_threads(nthreads);
  }
  for(sig = 0; sig < nsigma; sig++) {
    si = s[sig];
    Hhat = make_Hhat(dist_de, si);
    hhat = make_hhat(dist_nu, si);
#pragma omp parallel for num_threads(nthreads) if (parallel)
    for(l = 0; l < nlambda; l++) {
      alpha.slice(sig).col(l) = compute_alpha(Hhat, hhat, lambda(l));
    }
  }
  List out = List::create(
    Named("alpha") = alpha,
    Named("lambda") = lambda,
    Named("sigma") = s
  );
  return out;
}

