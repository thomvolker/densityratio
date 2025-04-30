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
arma::vec kmm_unconstrained_alpha(const arma::mat& Kdn,
                                  const arma::mat& Kdd,
                                  const arma::mat& Kd,
                                  const int nnu,
                                  const int nde) {
  arma::mat L1 = Kd.t() * Kdd * Kd;
  L1.diag() += 0.0001;
  arma::vec L2 = sum(Kd.t() * Kdn, 1);
  arma::vec alpha = arma::solve(L1, L2);
  return nde / nnu * alpha;
}

//[[Rcpp::export]]
arma::vec kmm_constrained_alpha(const arma::mat& Kdn,
                                const arma::mat& Kdd,
                                const arma::mat& Kd,
                                const int nnu,
                                const int nde,
                                const List settings) {

  Environment osqp = Environment::namespace_env("osqp");
  Function osqp_solve = osqp["solve_osqp"];

  double eps = 1 / sqrt(nde);
  arma::mat Amat  = join_cols(mean(Kd, 0), arma::eye(Kd.n_cols, Kd.n_cols));
  arma::mat lower = join_cols(arma::ones(1,1) * (1-eps),
                              arma::zeros(Kd.n_cols, 1));
  arma::mat upper = join_cols(arma::ones(1,1) * (1+eps),
                              100 * arma::ones(Kd.n_cols, 1));

  arma::mat P = Kd.t() * Kdd * Kd / (nde*nde);
  arma::mat q = -Kd.t() * sum(Kdn, 1) / (nde*nnu);

  List osqp_res = osqp_solve(P, q, Amat, lower, upper, settings);

  return osqp_res["x"];
}

// [[Rcpp::export]]
double kmm_cv_loss(const arma::mat& Kdn,
                   const arma::mat& Kdd,
                   const arma::mat& Kd,
                   const arma::mat& Kn,
                   const int& nfolds,
                   const arma::uvec& cv_ind_nu,
                   const arma::uvec& cv_ind_de,
                   bool constrained,
                   List settings) {

  double loss = 0;

  for (int k = 0; k < nfolds; k++) {

    arma::uvec idx_train_de = find(cv_ind_de != k);
    arma::uvec idx_train_nu = find(cv_ind_nu != k);
    arma::uvec idx_test_de = find(cv_ind_de == k);
    arma::uvec idx_test_nu = find(cv_ind_nu == k);
    arma::vec alpha;

    if (constrained) {
      alpha = kmm_constrained_alpha(Kdn.submat(idx_train_de, idx_train_nu),
                                    Kdd.submat(idx_train_de, idx_train_de),
                                    Kd.rows(idx_train_de),
                                    idx_train_nu.n_elem,
                                    idx_train_de.n_elem,
                                    settings);
    } else {
      alpha = kmm_unconstrained_alpha(Kdn.submat(idx_train_de, idx_train_nu),
                                Kdd.submat(idx_train_de, idx_train_de),
                                Kd.rows(idx_train_de),
                                idx_train_nu.n_elem,
                                idx_train_de.n_elem);
    }

    arma::vec r_nu = Kn.rows(idx_test_nu) * alpha;
    arma::vec r_de = Kd.rows(idx_test_de) * alpha;

    loss += mean(square(r_de) - 2 * mean(r_nu));
  }
  return loss / nfolds;
}


//[[Rcpp::export]]
List compute_kmm(const arma::mat& nu,
                 const arma::mat& de,
                 const arma::mat& ce,
                 const arma::mat& Dd,
                 const arma::vec& sigma,
                 const arma::uvec& cv_ind_nu,
                 const arma::uvec& cv_ind_de,
                 const bool& parallel,
                 int nthreads,
                 const bool& progressbar,
                 const bool& constrained,
                 const List settings) {

  double si;
  bool stopped = false;
  int nsigma = sigma.size();
  int nnu = nu.n_rows;
  int nde = de.n_rows;
  int nce = Dd.n_cols;
  int nfolds = max(cv_ind_nu) + 1;
  arma::mat Dn;

  arma::vec loss(nsigma);
  arma::mat Ddd = distance(de, de, false);
  arma::mat Ddn = distance(de, nu, false);
  arma::vec current_alpha(nce);
  arma::mat alpha(nce, nsigma);

  if (nfolds > 1) {
    Dn = distance(nu, ce, false);
  }

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
  Progress p(nsigma, progressbar);

  #ifdef _OPENMP
  #pragma omp parallel for num_threads(nthreads) private(si, current_alpha) if (parallel)
  #endif
    for (int sig = 0; sig < nsigma; sig++) {
      if(Progress::check_abort()) {
        stopped = true;
      } else {
        p.increment();
        si = sigma[sig];
        arma::mat Kdn = kernel_gaussian(Ddn, si);
        arma::mat Kdd = kernel_gaussian(Ddd, si);
        arma::mat Kd = kernel_gaussian(Dd, si);
        // constrained is always ran sequentially here; parallel = false in R code
        if (constrained) {
          current_alpha = kmm_constrained_alpha(Kdn, Kdd, Kd, nnu, nde, settings);
          alpha.col(sig) = current_alpha;
        } else {
          current_alpha = kmm_unconstrained_alpha(Kdn, Kdd, Kd, nnu, nde);
          alpha.col(sig) = current_alpha;
        }
        if (nfolds > 1) {
          arma::mat Kn = kernel_gaussian(Dn, si);
          loss[sig] = kmm_cv_loss(Kdn, Kdd, Kd, Kn, nfolds, cv_ind_nu, cv_ind_de, constrained, settings);
        }
      }
    }

  if (stopped) {
    if (progressbar) {
      Rprintf("\n");
    }
    Rcpp::stop("User terminated execution.");
    R_FlushConsole();
  }
  List out = List::create(
    Named("alpha") = alpha,
    Named("loss") = loss
  );
  return out;
}
