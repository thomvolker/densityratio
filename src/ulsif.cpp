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
arma::vec ulsif_compute_alpha(arma::mat Hhat, arma::vec hhat, double lambda, bool intercept) {
  if (intercept) {
    for (arma::uword i = 1; i < Hhat.n_cols; ++i) {
      Hhat(i, i) += lambda;
    }
  } else {
    Hhat.diag() += lambda;
  }
  arma::vec alpha = arma::solve(Hhat, hhat);
  return alpha;
}

//[[Rcpp::export]]
int set_threads(int nthreads) {

  #ifdef _OPENMP
  int max_threads = std::min(omp_get_max_threads(), omp_get_thread_limit());

  #else
  int max_threads = 1;
  #endif
  
  if (nthreads > max_threads) {
    nthreads = max_threads;
    std::string warn = "'nthreads' exceeds the maximum number of threads to use for parallel processing; it is set to the maximum number of threads available (" + std::to_string(nthreads) + ")";
    Rf_warningcall(R_NilValue, "%s", warn.c_str());
  }

  if (nthreads == 0) {
    nthreads = max_threads;
  }
  return nthreads;
}

//[[Rcpp::export]]
double compute_ulsif_loocv(const arma::mat& Gnu, const arma::mat& Gde, const arma::vec& hhat, double lambda, double a, bool intercept, const int& nnu, const int& nde, const arma::mat& Knu_nmin, const arma::mat& Kde_nmin) {
  
  double c1 = a / (nnu - 1.0);
  double c2 = (1.0 - a) / (nde - 1.0);
  
  arma::mat A = c1 * Gnu + c2 * Gde;
  A.diag() += lambda;
  if (intercept) {
    A(0, 0) -= lambda;
  }

  arma::mat Ainv = arma::inv(A);
  arma::mat P = Knu_nmin * Ainv;
  arma::mat Q = Kde_nmin * Ainv;

  arma::vec s11 = arma::sum(P % Knu_nmin, 1);
  arma::vec s12 = arma::sum(P % Kde_nmin, 1);
  arma::vec s22 = arma::sum(Q % Kde_nmin, 1);

  arma::vec t1 = (nnu * (P * hhat) - s11) / (nnu - 1.0);
  arma::vec t2 = (nnu * (Q * hhat) - s12) / (nnu - 1.0);

  arma::vec d11 = 1.0 - c1 * s11;
  arma::vec d22 = 1.0 - c2 * s22;
  arma::vec det = d11 % d22 - c1 * c2 * arma::square(s12);

  arma::vec wnu = (d22 % t1 + c2 * (s12 % t2)) / det;
  arma::vec wde = (c1 * (s12 % t1) + d11 % t2) / det;

  return arma::mean(0.5 * a * arma::square(wnu) + 0.5 * (1.0 - a) * arma::square(wde) - wnu);
}

//[[Rcpp::export]]
List compute_ulsif(const arma::mat& dist_nu, const arma::mat& dist_de, const bool& intercept, const arma::vec& sigma, const arma::vec& lambda, const bool& parallel, int nthreads, const bool& progressbar, const double& a = 0) {

  double si, la;
  bool stopped = false;
  int nsigma  = sigma.size();
  int nlambda = lambda.size();
  int nnu = dist_nu.n_rows;
  int nde = dist_de.n_rows;
  int nmin = (nnu <= nde) ? nnu : nde;
  int ncol = dist_nu.n_cols;
  int sig;

  arma::mat loocv(nsigma, nlambda);
  arma::mat Hhat(ncol, ncol);
  arma::mat Gnu(ncol, ncol, arma::fill::zeros);
  arma::mat Gde(ncol, ncol);
  arma::vec hhat(ncol);
  arma::mat Knu(nnu, ncol);
  arma::mat Kde(nde, ncol);
  arma::vec current_alpha(ncol);
  arma::cube alpha(ncol, nsigma, nlambda);

  #ifdef _OPENMP
  if (parallel) {
    nthreads = set_threads(nthreads);
  }
  #else
  if (parallel) {
    std::string warn = "OpenMP is not available, parallel processing is disabled.";
    Rf_warningcall(R_NilValue, "%s", warn.c_str());
  }
  #endif
  Progress p(nsigma * nlambda, progressbar);

  for(sig = 0; sig < nsigma; sig++) {
    si = sigma[sig];
    Knu = kernel_gaussian(dist_nu, si);
    Kde = kernel_gaussian(dist_de, si);
    
    Gde = Kde.t() * Kde;
    Hhat = (1-a) * Gde / nde;
    if (a > 0) {
      Gnu = Knu.t() * Knu;
      Hhat += a * Gnu / nnu;
    }

    hhat = arma::mean(Knu, 0).t();
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(nthreads) private(la, current_alpha) if (parallel)
  #endif
    for(int l = 0; l < nlambda; l++) {
      if(Progress::check_abort()) {
        #pragma omp critical
        { 
          stopped = true; 
        }
      } else {
        la = lambda[l];
        p.increment();
        current_alpha = ulsif_compute_alpha(Hhat, hhat, la, intercept);
        alpha.slice(l).col(sig) = current_alpha;
        loocv(sig, l) = compute_ulsif_loocv(Gnu, Gde, hhat, la, a, intercept, nnu, nde, Knu.rows(0, nmin-1), Kde.rows(0,nmin-1));
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

