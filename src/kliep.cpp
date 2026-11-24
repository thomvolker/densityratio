//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>
//#include <omp.h>
#include "densityratio.h"
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::mat make_Phi(const arma::mat& dist_nu, double sigma) {
  arma::mat Phi = kernel_gaussian(dist_nu, sigma);
  return Phi;
}

//[[Rcpp::export]]
arma::mat make_phibar(const arma::mat& dist_de, double sigma) {
  arma::mat phibar = kernel_gaussian(dist_de, sigma);
  return mean(phibar, 0).t();
}

arma::vec kliep_compute_alpha(const arma::mat& Phi, const arma::vec& phibar, double tol, int maxit, bool progressbar) {

  bool stopped = false;
  bool conv = false;
  int iter = 0;
  double s_new, score;

  arma::vec alpha = arma::ones<arma::vec>(Phi.n_cols);
  arma::vec alpha_new = arma::ones<arma::vec>(Phi.n_cols);
  arma::vec ratio = arma::zeros<arma::vec>(Phi.n_cols);
  arma::vec Phi_alpha = arma::zeros<arma::vec>(Phi.n_rows);
  arma::vec Phi_alpha_new = arma::zeros<arma::vec>(Phi.n_rows);

  alpha /= dot(phibar, alpha);
  Phi_alpha = clamp(Phi * alpha, 1e-300, datum::inf);
  score = mean(log(Phi_alpha));

  while(!conv) {
    iter++;
    if (iter % 100 == 0 && Progress::check_abort()) {
      stopped = true;
      break;
    }
    // gradient of mean(log(Phi * alpha)), divided by phibar
    ratio = (Phi.t() * (1 / Phi_alpha)) / (Phi.n_rows * phibar);
    // multiplicative (EM) step: keeps alpha >= 0 and phibar' alpha = 1
    alpha_new = alpha % ratio;
    alpha_new /= dot(phibar, alpha_new);
    Phi_alpha_new = clamp(Phi * alpha_new, 1e-300, datum::inf);
    s_new = mean(log(Phi_alpha_new));
    if (!(s_new > score)) {
      conv = true;
    } else {
      score = s_new;
      alpha = alpha_new;
      Phi_alpha = Phi_alpha_new;
      // log(max(ratio)) bounds the distance to the optimal objective value
      conv = log(ratio.max()) < tol || iter == maxit;
    }
  }
  if (stopped) {
    if (progressbar) {
      Rprintf("\n");
    }
    Rcpp::stop("User terminated execution.");
    R_FlushConsole();
  }
  return alpha;
}

//[[Rcpp::export]]
List compute_kliep(const arma::mat& dist_nu, const arma::mat& dist_de, const arma::vec& sigma, const double& tol, const int& maxit, arma::vec cv_ind, bool progressbar) {

  int dim = dist_nu.n_cols;
  int nsigma = sigma.size();
  double sig;
  int nfold = max(cv_ind) + 1;
  Progress p(nsigma*nfold, progressbar);

  arma::mat alpha = arma::ones<arma::mat>(dim, nsigma);
  arma::mat Phi = arma::zeros<arma::mat>(dim, dim);
  arma::vec phibar = arma::zeros<arma::vec>(dim);
  arma::vec alpha_tmp;
  arma::vec cv_score = arma::zeros<arma::vec>(nsigma);


  for (int s = 0; s < nsigma; s++) {
    sig = sigma[s];
    Phi = make_Phi(dist_nu, sig);
    phibar = make_phibar(dist_de, sig);

    if (nfold > 1) {
      for (int i = 0; i < nfold; i++) {
        p.increment();
        alpha_tmp = kliep_compute_alpha(Phi.rows(find(cv_ind != i)), phibar, tol, maxit, progressbar);
        cv_score(s) += mean(log(Phi.rows(find(cv_ind == i)) * alpha_tmp));
      }
    } else {
      p.increment();
    }
    alpha.col(s) = kliep_compute_alpha(Phi, phibar, tol, maxit, progressbar);
  }
  List out = List::create(
    Named("alpha") = alpha,
    Named("cv_score") = cv_score
  );
  return out;
}