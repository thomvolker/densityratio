#ifndef DENSITYRATIO
#define DENSITYRATIO

#include <Rcpp.h>

arma::mat distance(arma::mat X, arma::mat Y, bool symmetric = false);
arma::mat kernel_gaussian(arma::mat dist, double sigma);
Rcpp::List compute_ulsif(arma::mat dist_nu, arma::mat dist_de, arma::vec sigma, arma::vec lambda, bool parallel, int nthreads, bool progressbar);

#endif

