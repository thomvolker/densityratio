#ifndef DENSITYRATIO
#define DENSITYRATIO

#include <Rcpp.h>
#include <RcppArmadillo.h>

arma::mat  distance(const arma::mat& X, const arma::mat& Y, const bool& intercept);
arma::mat  kernel_gaussian(const arma::mat& dist, double sigma);
Rcpp::List compute_ulsif(const arma::mat& dist_nu, const arma::mat& dist_de, const arma::vec& sigma, const arma::vec& lambda, const bool& parallel, const int& nthreads, const bool& progressbar);
double     compute_ulsif_loocv(arma::mat Hhat, arma::mat hhat, double lambda, const int& nnu, const int& nde, const int& nmin, const int& ncol, arma::mat Knu_nmin, arma::mat Kde_nmin);
int        set_threads(int nthreads);
arma::vec  ulsif_compute_alpha(arma::mat Hhat, arma::vec hhat, double lambda);

#endif

