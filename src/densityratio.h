#include <Rcpp.h>
#ifndef DENSITYRATIO
#define DENSITYRATIO

arma::mat distance(arma::mat X, arma::mat Y, bool symmetric = false);
arma::vec make_sigma(arma::mat dist, arma::vec sigma);
arma::mat kernel_gaussian(arma::mat dist, double sigma);

#endif

