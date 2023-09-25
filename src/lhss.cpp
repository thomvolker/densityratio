//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>
#include "densityratio.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat make_UV(arma::mat U) {
  arma::mat Q, R;
  qr(Q, R, U);
  return Q * sign(as_scalar(R(0,0)));
}

//[[Rcpp::export]]
arma::vec check_sigma_cpp(int nsigma, Nullable<arma::vec> sigma_quantile, Nullable<arma::vec> sigma, arma::mat dist) {
  Function check_sigma("check.sigma");
  NumericVector s = check_sigma(nsigma, sigma_quantile, sigma, dist);
  return as<arma::vec>(wrap(s));
}

// // [[Rcpp::export]]
// List compute_alpha_lhss(arma::mat nu, arma::mat de, arma::mat ce, int m,
//                         int nsigma, Nullable<arma::vec> sigma_quantile,
//                         Nullable<arma::vec> sigma, arma::vec lambda,
//                         int maxit, bool parallel, int nthreads, bool progressbar) {
//
//   Function check_sigma("check.sigma");
//
//   int ncol = nu.n_cols;
//   int n_nu = nu.n_rows;
//   int n_de = de.n_rows;
//   int n_ce = ce.n_rows;
//
//   arma::mat U = make_UV(arma::ones<arma::mat>(ncol, m)).cols(0, m-1);
//
//
//   nu = nu * U;
//   de = de * U;
//   ce = ce * U;
//
//   arma::mat dist_nu = distance(nu, ce);
//   arma::mat dist_de = distance(de, ce);
//
//   arma::vec s = check_sigma_cpp(nsigma, sigma_quantile, sigma, dist_nu);
//
//   int nlambda = lambda.size();
//   nsigma  = s.size();
//
//   List ulsif_init = compute_ulsif(dist_nu, dist_de, s, lambda, parallel, nthreads, false);
//   arma::vec loocv_scores = ulsif_init["loocv_score"];
//   arma::vec PD_scores = ulsif_init["PD"];
//   int min_score_index = loocv_scores.index_min();
//   int lambda_min_index = min_score_index % nlambda;
//   int sigma_min_index  = floor(min_score_index / nlambda);
//   arma::cube alpha = ulsif_init["alpha"];
//   arma::vec alpha_min = alpha.slice(sigma_min_index).col(lambda_min_index);
//   arma::vec alphah = alpha_min;
//   alphah.elem(find(alpha_min < 0)).zeros();
//
//   double sigma_opt = s(sigma_min_index);
//   double lambda_opt = lambda(lambda_min_index);
//   double PD_opt = as_scalar(PD_scores(min_score_index));
//   arma::mat U_opt = U;
//   double score_opt = loocv_scores(min_score_index);
//
//   int decPDcount = 0;
//   bool conv = false;
//   int iter = 0;
//
//   arma::vec one_n_nu = arma::ones<arma::vec>(n_nu);
//   arma::vec one_n_de = arma::ones<arma::vec>(n_de);
//   arma::vec one_n_ce_T = arma::ones<arma::vec>(n_ce).t();
//   arma::vec one_m_T = arma::ones<arma::vec>(m).t();
//   arma::mat UV = arma::zeros<arma::mat>(ncol, ncol);
//   arma::mat V = arma::zeros<arma::mat>(ncol, ncol - m);
//   arma::mat dist_nu_u = arma::zeros<arma::mat>(n_nu, n_ce);
//   arma::mat Knu = arma::zeros<arma::mat>(n_nu, n_ce);
//
//   while (!conv) {
//     iter++;
//
//     UV = make_UV(U);
//     U  = UV.cols(0, m-1);
//     V  = UV.cols(m, ncol - 1);
//
//     nu = nu * U;
//     de = de * U;
//     ce = ce * U;
//
//     dist_nu_u = distance(nu, ce);
//     s = check_sigma_cpp(nsigma, sigma_quantile, sigma, dist_nu);
//     Knu = kernel_gaussian(dist_nu_u, s);
//   }
//
//   List out = List::create(
//     Named("alpha") = alpha_min
//   );
//   return out;
// }
