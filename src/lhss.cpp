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
arma::mat make_UV(arma::mat U) {
  arma::mat Q, R;
  arma::qr(Q, R, U);
  return Q * sign(as_scalar(R(0,0)));
}

//[[Rcpp::export]]
double get_sigma_lhss(arma::mat dist, arma::vec sigma, bool quantiles) {
  double s;
  if (quantiles) {
    arma::vec nonzero_dist = nonzeros(dist);
    arma::vec temp_s = quantile(nonzero_dist, sigma);
    s = as_scalar(temp_s);
  } else {
    s = as_scalar(sigma);
  }
  return s;
}


// [[Rcpp::export]]
List lhss_compute_alpha(arma::mat nu, arma::mat de,
                        arma::mat ce, bool symmetric,
                        int m, bool intercept,
                        arma::vec sigma, bool quantiles,
                        arma::vec lambda, int maxit,
                        bool progressbar) {

  // initialize objects
  int p = nu.n_cols;
  bool is_subspace = m < p;
  int n_nu = nu.n_rows;
  int n_de = de.n_rows;
  int n_ce = ce.n_rows;
  int nbasis = intercept ? ce.n_rows + 1 : ce.n_rows; // check for intercept (for dimensions later)
  int nlambda = lambda.size();
  int nsigma  = sigma.size();
  int nmin = (n_nu <= n_de) ? n_nu : n_de;
  double s, PD, PD_opt;
  bool conv;
  bool stopped = false;
  int PDcount, iter;

  // initialize empty matrices
  arma::mat UV(p, p);
  arma::mat eyemat = arma::eye(p, m);
  arma::mat nu_u(n_nu, m);
  arma::mat de_u(n_de, m);
  arma::mat ce_u(n_ce, m);
  arma::mat dist_nu_u(n_nu, nbasis);
  arma::mat dist_de_u(n_de, nbasis);
  arma::mat Knu(n_nu, nbasis);
  arma::mat Kde(n_de, nbasis);
  arma::mat Knu_nmin(nmin, nbasis);
  arma::mat Kde_nmin(nmin, nbasis);
  arma::mat Hhat(nbasis, nbasis);
  arma::mat Hhat_opt(nbasis, nbasis);
  arma::mat dPd1(p, m);
  arma::mat dPd2(p, m);
  arma::mat dPd(m, p);
  arma::mat dM(p, p);
  arma::mat eM(p, p);
  arma::mat temp11(n_nu, m);
  arma::mat temp12(n_nu, p);
  arma::mat temp21(n_de, m);
  arma::mat temp22(n_de, p);
  arma::vec Hhat_diag(nbasis);
  arma::vec hhat(nbasis);
  arma::vec hhat_opt(nbasis);
  arma::vec current_alpha(nbasis);
  arma::vec current_sigma(1);
  arma::cube alpha(nbasis, nsigma, nlambda);
  arma::mat sigmaopt(nsigma, nlambda);
  arma::mat loocv(nsigma, nlambda);
  arma::cube Umat(p, m, nsigma * nlambda);

  // initialize progressbar
  Progress prog(nsigma * nlambda, progressbar);

  // start main outer loop over sigma
  for (int sig = 0; sig < nsigma; sig++) {
    current_sigma = sigma[sig];                // obtain current sigma (can be a quantile)
    // start inner loop over lambda
    for (int lam = 0; lam < nlambda; lam++) {

      if (Progress::check_abort()) { // check whether program is stopped
        stopped = true;
      } else {
        prog.increment();

        UV = make_UV(arma::ones<arma::mat>(p, m)); // initialize subspace
        nu_u = nu * UV.cols(0, m-1);
        de_u = de * UV.cols(0, m-1);
        ce_u = ce * UV.cols(0, m-1);

        dist_nu_u = distance(nu_u, ce_u, intercept); // calculate distance matrices
        dist_de_u = distance(de_u, ce_u, intercept);

        s = get_sigma_lhss(dist_nu_u, current_sigma, quantiles); // extract actual sigma from distances if sigma is a quantile

        Knu = kernel_gaussian(dist_nu_u, s); // calculate kernel matrices
        Kde = kernel_gaussian(dist_de_u, s);
        Hhat = Kde.t() * Kde / n_de;
        hhat = arma::mean(Knu, 0).t();

        // compute starting lambda
        current_alpha = ulsif_compute_alpha(Hhat, hhat, lambda(lam));
        // set negative elements to zero
        current_alpha.elem(find(current_alpha < 0)).zeros();
        // calculate pearson divergence
        PD_opt = dot(hhat, current_alpha) - 1/2;

        Umat.slice(lam * nsigma + sig) = UV.cols(0, m-1);
        alpha.slice(lam).col(sig) = current_alpha;
        sigmaopt(sig, lam) = s;
        Hhat_opt = Hhat;
        hhat_opt = hhat;
        Knu_nmin = Knu.rows(0, nmin - 1);
        Kde_nmin = Kde.rows(0, nmin - 1);

        PDcount = 0;
        conv    = false;
        iter    = 0;

        // optimize subspace
        while (!conv) {
          iter++;

          // update U given alpha {
          UV = make_UV(UV.cols(0, m-1));

          nu_u = nu * UV.cols(0, m-1);
          de_u = de * UV.cols(0, m-1);
          ce_u = ce * UV.cols(0, m-1);

          dPd1.zeros();
          dPd2.zeros();

          dist_nu_u = distance(nu_u, ce_u, intercept);
          dist_de_u = distance(de_u, ce_u, intercept);

          s = get_sigma_lhss(dist_nu_u, current_sigma, quantiles);

          Knu = kernel_gaussian(dist_nu_u, s);
          Kde = kernel_gaussian(dist_de_u, s);

          for (int i = intercept ? 1 : 0; i < n_ce; i++) {
            temp11 = -(nu_u - repmat(ce_u.row(i), n_nu, 1)) % repmat(Knu.col(i), 1, m);
            temp12 = nu - repmat(ce.row(i), n_nu, 1);
            dPd1 += temp12.t() * temp11 * current_alpha(i);

            temp21 = (de_u - repmat(ce_u.row(i), n_de, 1)) % repmat(Kde.col(i), 1, m);
            temp22 = de - repmat(ce.row(i), n_de, 1);

            for (int j = 0; j < m; j++) {
              dPd2.col(j) -= temp22.t() * (repmat(temp21.col(j), 1, nbasis) % Kde) * current_alpha * current_alpha(i) * 2;
            }
          }
          dPd = (dPd1 / n_nu / s - dPd2 / n_de / s / 2).t();
          if (is_subspace) {
            dM.submat(0, m, m-1, p - 1) = dPd * UV.cols(m, p-1);
            dM.submat(m, 0, p-1, m-1) = - (dPd * UV.cols(m, p-1)).t();
          }
          eM = expmat(dM / maxit * (maxit - iter + 1));
          UV.cols(0, m-1) = UV * eM.t() * eyemat;
          // } end update U given alpha
          // update alpha given U {
          nu_u = nu * UV.cols(0, m-1);
          de_u = de * UV.cols(0, m-1);
          ce_u = ce * UV.cols(0, m-1);

          dist_nu_u = distance(nu_u, ce_u, intercept);
          dist_de_u = distance(de_u, ce_u, intercept);

          Knu = kernel_gaussian(dist_nu_u, s);
          Kde = kernel_gaussian(dist_de_u, s);
          Hhat = Kde.t() * Kde / n_de;
          hhat = mean(Knu, 0).t();


          current_alpha = ulsif_compute_alpha(Hhat, hhat, lambda(lam));
          current_alpha.elem(find(current_alpha < 0)).zeros();
          // } end update alpha given U
          PD = dot(hhat, current_alpha) - 1/2;

          // check whether we're heading in the right direction
          if (PD <= PD_opt) {
            PDcount++;
          } else {
            PD_opt = PD;
            PDcount = PDcount - 1 < 0 ? 0 : PDcount - 1;
            Umat.slice(lam * nsigma + sig) = UV.cols(0, m-1);
            alpha.slice(lam).col(sig) = current_alpha;
            sigmaopt(sig, lam) = s;
            Hhat_opt = Hhat;
            hhat_opt = hhat;
            Knu_nmin = Knu.rows(0, nmin - 1);
            Kde_nmin = Kde.rows(0, nmin - 1);
          }
          if ((iter == maxit) | (PDcount > 20)) {
            conv = true;
          }
        }
        loocv(sig, lam) = compute_ulsif_loocv(Hhat_opt, hhat_opt, lambda(lam), n_nu, n_de, nmin, nbasis, Knu_nmin, Kde_nmin);
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
  List out = List::create(
    Named("alpha") = alpha,
    Named("Uopt") = Umat,
    Named("sigmaopt") = sigmaopt,
    Named("loocv") = loocv
  );
  return out;
}
