//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>
#include "densityratio.h"
#include <progress.hpp>
#include <progress_bar.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;


//[[Rcpp::export]]
arma::vec ulsif_compute_alpha(arma::mat Hhat, const arma::vec& hhat, const double& lambda) {
  Hhat.diag() += lambda;
  arma::vec alpha = arma::solve(Hhat, hhat);
  return alpha;
}

//[[Rcpp::export]]
int set_threads(int nthreads) {
  #ifdef _OPENMP
  int max_threads = omp_get_max_threads();
  #else
  int max_threads = 1;
  #endif
  if (nthreads == 0) {
    nthreads = max_threads;
  }
  if (nthreads > max_threads) {
    nthreads = max_threads;
    std::string warn = "'nthreads' exceeds the maximum number of threads to use for parallel processing; it is set to the maximum number of threads available (" + std::to_string(nthreads) + ")";
    Rf_warningcall(R_NilValue, warn.c_str());
  }
  return nthreads;
}

//[[Rcpp::export]]
double compute_ulsif_loocv(arma::mat Hhat, const arma::mat& hhat, const double& lambda, const int& nnu, const int& nde, const int& nmin, const int& ncol, const arma::mat& Knu_nmin, const arma::mat& Kde_nmin) {
  double la = lambda * (nde - 1) / nde;

  arma::vec onen = arma::ones<arma::vec>(nmin);
  arma::vec oneb = arma::ones<arma::vec>(ncol);

  arma::mat Binv = arma::inv(Hhat + arma::diagmat(oneb * la));
  arma::mat B0   = repmat(Binv * hhat, 1, nmin);
  arma::mat B02  = Kde_nmin * Binv;
  B0 += B02.t() * diagmat(B02 * hhat / (nde - sum(Kde_nmin % B02, 1)));

  arma::mat B1  = (Knu_nmin * Binv).t();
  arma::vec B12 = sum(Knu_nmin % B02, 1) / (nmin - sum(Kde_nmin % B02, 1));
  B1 += B02.t() * diagmat(B12);

  arma::mat B2  = (nnu * B0 - B1) * (nde - 1) / (nde * (nnu - 1));
  B2.elem(find(B2 < 0)).zeros();

  arma::vec wtr = sum(Kde_nmin % B2.t(), 1);
  arma::vec wte = sum(Knu_nmin % B2.t(), 1);

  return as_scalar(arma::dot(wtr, wtr) / (2*nmin) - mean(wte));
}

//[[Rcpp::export]]
List compute_ulsif(arma::mat dist_nu, arma::mat dist_de, arma::vec sigma, arma::vec lambda, bool parallel, int nthreads, bool progressbar) {

  double si, la;
  bool stopped = false;
  int nsigma  = sigma.size();
  int nlambda = lambda.size();
  int nnu = dist_nu.n_rows;
  int nde = dist_de.n_rows;
  int nmin = (nnu <= nde) ? nnu : nde;
  int ncol = dist_nu.n_cols;
  int sig;

  arma::mat loocv = arma::zeros<arma::mat>(nsigma, nlambda);
  arma::mat Hhat = arma::zeros<arma::mat>(ncol, ncol);
  arma::vec hhat = arma::zeros<arma::vec>(ncol);
  arma::mat Knu = arma::zeros<arma::mat>(nnu, ncol);
  arma::mat Kde = arma::zeros<arma::mat>(nde, ncol);
  arma::vec current_alpha = arma::zeros<arma::vec>(ncol);
  arma::cube alpha(ncol, nsigma, nlambda, arma::fill::zeros);

  #ifdef _OPENMP
  if (parallel) {
    nthreads = set_threads(nthreads);
  }
  #else
  if (parallel) {
    std::string warn = "OpenMP is not available, parallel processing is disabled.";
    Rf_warningcall(R_NilValue, warn.c_str());
  }
  #endif
  Progress p(nsigma * nlambda, progressbar);

  for(sig = 0; sig < nsigma; sig++) {
    si = sigma[sig];
    Knu = kernel_gaussian(dist_nu, si);
    Kde = kernel_gaussian(dist_de, si);
    Hhat = Kde.t() * Kde / nde;
    hhat = arma::mean(Knu, 0).t();
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(nthreads) private(la, current_alpha) if (parallel)
  #endif
    for(int l = 0; l < nlambda; l++) {
      if(Progress::check_abort()) {
        stopped = true;
      } else {
        la = lambda[l];
        p.increment();
        current_alpha =  ulsif_compute_alpha(Hhat, hhat, la);
        alpha.slice(l).col(sig) = current_alpha;
        loocv(sig, l) = compute_ulsif_loocv(Hhat, hhat, la, nnu, nde, nmin, ncol, Knu.rows(0, nmin-1), Kde.rows(0,nmin-1));
      }
    }
    if (stopped) {
      if (progressbar) {
        Rprintf("\n");
      }
      Rcpp::stop("User terminated execution.");
      R_FlushConsole();
    }
  }

  List out = List::create(
    Named("alpha") = alpha,
    Named("loocv_score") = loocv
  );
  return out;
}
