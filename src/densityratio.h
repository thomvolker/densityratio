#ifndef DENSITYRATIO
#define DENSITYRATIO

#include <Rcpp.h>
#include <RcppArmadillo.h>

arma::mat  distance(const arma::mat& X, const arma::mat& Y, const bool& intercept);
arma::mat  kernel_gaussian(const arma::mat& dist, double sigma);
Rcpp::List compute_ulsif(const arma::mat& dist_nu, const arma::mat& dist_de, const bool& intercept, const arma::vec& sigma, const arma::vec& lambda, const bool& parallel, int nthreads, const bool& progressbar, const double& a);
double     compute_ulsif_loocv(const arma::mat& Gnu, const arma::mat& Gde, const arma::vec& hhat, double lambda, double a, bool intercept, const int& nnu, const int& nde, const arma::mat& Knu_nmin, const arma::mat& Kde_nmin);
int        set_threads(int nthreads);
arma::vec  ulsif_compute_alpha(arma::mat Hhat, arma::vec hhat, double lambda, bool intercept);

#endif

