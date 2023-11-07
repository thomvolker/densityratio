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
  arma::mat Phi = kernel_gaussian(dist_nu, sigma).t();
  return Phi;
}

//[[Rcpp::export]]
arma::mat make_phibar(const arma::mat& dist_de, double sigma) {
  arma::mat phibar = kernel_gaussian(dist_de, sigma);
  return mean(phibar, 0).t();
}

//[[Rcpp::export]]
arma::vec kliep_compute_alpha(const arma::mat& Phi, const arma::vec& phibar, const arma::vec& phibar_corr, const arma::vec& epsilon, int nepsilon, int maxit, bool progressbar) {

  bool stopped = false;
  bool conv = false;
  int iter = 0;
  double s_new, e, score;

  arma::vec alpha = arma::ones<arma::vec>(Phi.n_rows);
  arma::vec alpha_new = arma::ones<arma::vec>(Phi.n_rows);
  arma::vec Phi_t_alpha = arma::zeros<arma::vec>(Phi.n_cols);
  arma::vec Phi_t_alpha_new = arma::zeros<arma::vec>(Phi.n_cols);

  alpha += phibar_corr * as_scalar(1 - arma::dot(phibar, alpha));
  alpha.elem(find(alpha < 0)).zeros();
  alpha /= dot(phibar, phibar);
  score = mean(log(Phi.t() * alpha));

  for(int eps = 0; eps < nepsilon; eps++) {
    conv = false;
    iter = 0;
    e = epsilon[eps];

    if(Progress::check_abort()) {
      stopped = true;
    } else {
      while(!conv) {
        iter++;
        Phi_t_alpha = Phi.t() * alpha;
        Phi_t_alpha.replace(0, sqrt(datum::eps));
        alpha_new = alpha + (e * Phi) * (1 / (Phi_t_alpha));
        alpha_new += phibar_corr * as_scalar(1 - arma::dot(phibar, alpha_new));
        alpha_new.elem(find(alpha_new < 0)).zeros();
        alpha_new /= dot(phibar, alpha_new);
        Phi_t_alpha_new = Phi.t() * alpha_new;
        Phi_t_alpha_new.replace(0, sqrt(datum::eps));
        s_new = mean(log((Phi_t_alpha_new)));
        if (s_new <= score || iter == maxit) {
          conv = true;
        } else {
          score = s_new;
          alpha = alpha_new;
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
  }
  return alpha;
}

//[[Rcpp::export]]
List compute_kliep(const arma::mat& dist_nu, const arma::mat& dist_de, const arma::vec& sigma, const arma::vec& epsilon, const int& maxit, arma::vec cv_ind, bool progressbar) {

  int dim = dist_nu.n_cols;
  int nsigma = sigma.size();
  int nepsilon = epsilon.size();
  double sig;
  int nfold = max(cv_ind) + 1;
  Progress p(nsigma*nfold, progressbar);

  arma::mat alpha = arma::ones<arma::mat>(dim, nsigma);
  arma::mat Phi = arma::zeros<arma::mat>(dim, dim);
  arma::vec phibar = arma::zeros<arma::vec>(dim);
  arma::vec phibar_corr = arma::zeros<arma::vec>(dim);
  arma::vec alpha_tmp;
  arma::vec cv_score = arma::zeros<arma::vec>(nsigma);


  for (int s = 0; s < nsigma; s++) {
    sig = sigma[s];
    Phi = make_Phi(dist_nu, sig);
    phibar = make_phibar(dist_de, sig);
    phibar_corr = phibar / arma::dot(phibar, phibar);

    if (nfold > 1) {
      for (int i = 0; i < nfold; i++) {
        p.increment();
        alpha_tmp = kliep_compute_alpha(Phi.cols(find(cv_ind != i)), phibar, phibar_corr, epsilon, nepsilon, maxit, progressbar);
        cv_score(s) += mean(log(Phi.cols(find(cv_ind == i)).t() * alpha_tmp));
      }
    } else {
      p.increment();
    }
    alpha.col(s) = kliep_compute_alpha(Phi, phibar, phibar_corr, epsilon, nepsilon, maxit, progressbar);
  }
  List out = List::create(
    Named("alpha") = alpha,
    Named("cv_score") = cv_score
  );
  return out;
}


