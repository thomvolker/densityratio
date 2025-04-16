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

// [[Rcpp::export]]
arma::mat compute_psihat(const arma::mat& K,
                         const arma::mat& Evecs,
                         const arma::vec& Evals,
                         const int& maxJ,
                         const int& ncol) {
  return K * Evecs.cols(ncol - maxJ, ncol - 1) * diagmat(sqrt(ncol)/Evals.subvec(ncol - maxJ, ncol - 1));
}

// [[Rcpp::export]]
arma::rowvec spectral_cv_loss(const arma::mat& Knu,
                           const arma::mat& Kde,
                           const arma::vec& m,
                           const int& maxM,
                           const int& nfolds,
                           const arma::uvec& cv_ind_nu,
                           const arma::uvec& cv_ind_de,
                           const int& nthreads,
                           const bool& parallel) {

  arma::vec loss(m.size());

  #ifdef _OPENMP
  #pragma omp parallel for num_threads(nthreads) if (parallel)
  #endif
  for (int k = 0; k < nfolds; k++) {

    arma::uvec idx_train_de = find(cv_ind_de != k);
    arma::uvec idx_train_nu = find(cv_ind_nu != k);
    arma::uvec idx_test_de = find(cv_ind_de == k);
    arma::uvec idx_test_nu = find(cv_ind_nu == k);

    arma::vec EigVals;
    arma::mat EigVecs;

    eig_sym(EigVals, EigVecs, Kde.submat(idx_train_de, idx_train_de));

    int ncol = EigVecs.n_cols;

    arma::vec betatilde = mean(compute_psihat(Knu.submat(idx_train_nu, idx_train_de), EigVecs, EigVals, maxM, ncol), 0).t();

    arma::mat psitest_nu = compute_psihat(Knu.submat(idx_test_nu, idx_train_de), EigVecs, EigVals, maxM, ncol);
    arma::mat psitest_de = compute_psihat(Kde.submat(idx_test_de, idx_train_de), EigVecs, EigVals, maxM, ncol);

    for (arma::uword j = 0; j < m.size(); j++) {
      int start = maxM - m(j);
      arma::vec beta_nu = psitest_nu.cols(start, maxM-1) * betatilde.subvec(start, maxM-1);
      arma::vec beta_de = psitest_de.cols(start, maxM-1) * betatilde.subvec(start, maxM-1);

      beta_nu.elem(find(beta_nu < 0)).zeros();
      beta_de.elem(find(beta_de < 0)).zeros();

      loss(j) += as_scalar(mean(square(beta_de)) - 2 * mean(beta_nu));
    }
  }
  return loss.t() / nfolds;
}



// [[Rcpp::export]]
List spectral_dre(const arma::mat& dist_nu,
                  const arma::mat& dist_de,
                  const arma::vec& m,
                  const arma::vec& sigma,
                  const arma::uvec& cv_ind_nu,
                  const arma::uvec& cv_ind_de,
                  const bool& parallel,
                  int nthreads,
                  const bool& progressbar) {

  int ncol = dist_de.n_cols;
  int nsigma = sigma.size();
  int nsubspace = m.size();
  int maxM = max(m);
  int nfolds = max(cv_ind_nu) + 1;
  arma::mat loss(nsigma, nsubspace);
  arma::mat betatilde(maxM, nsigma);
  arma::mat Evals(ncol, nsigma);
  arma::cube Evecs(ncol, ncol, nsigma);

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

  for (int sig = 0; sig < nsigma; sig++) {
    if(Progress::check_abort()) {
      if (progressbar) Rprintf("\n");
      Rcpp::stop("User terminated execution.");
      R_FlushConsole();
    } else {
      p.increment();
      double si = sigma(sig);
      arma::mat Knu = kernel_gaussian(dist_nu, si);
      arma::mat Kde = kernel_gaussian(dist_de, si);

      arma::vec Evals_sig;
      arma::mat Evecs_sig;
      eig_sym(Evals_sig, Evecs_sig, Kde);

      Evecs.slice(sig) = Evecs_sig;
      Evals.col(sig) = Evals_sig;

      arma::mat psihat = compute_psihat(Knu, Evecs.slice(sig), Evals.col(sig), maxM, ncol);
      betatilde.col(sig) = mean(psihat, 0).t();

      if (nfolds > 1) {
        loss.row(sig) = spectral_cv_loss(Knu, Kde, m, maxM, nfolds, cv_ind_nu, cv_ind_de, nthreads, parallel);
      }
    }
  }

  return List::create(Named("betatilde") = betatilde,
                      Named("loss") = loss,
                      Named("Evecs") = Evecs,
                      Named("Evals") = Evals);
}
