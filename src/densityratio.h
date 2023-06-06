#include <Rcpp.h>
#ifndef DENSITYRATIO
#define DENSITYRATIO

arma::mat distance(arma::mat X, arma::mat Y, bool symmetric = false);
arma::mat kernel_gaussian(arma::mat dist, double sigma);

#endif

