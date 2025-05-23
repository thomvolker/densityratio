# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Create a Gram matrix with squared Euclidean distances between
#' observations in the input matrix \code{X} and the input matrix \code{Y}
#' @name distance
#' @param X A numeric input matrix
#' @param Y A numeric input matrix with the same variables as \code{X}
#' @param intercept Logical indicating whether an intercept should be added to
#' the estimation procedure. In this case, the first column is an all-zero
#' column (which will be transformed into an all-ones column in the kernel).
NULL

#' Create gaussian kernel gram matrix from distance matrix
#' @name kernel_gaussian
#' @param dist A numeric distance matrix
#' @param sigma A scalar with the length-scale parameter
NULL

distance <- function(X, Y, intercept = FALSE) {
    .Call(`_densityratio_distance`, X, Y, intercept)
}

kernel_gaussian <- function(dist, sigma) {
    .Call(`_densityratio_kernel_gaussian`, dist, sigma)
}

make_Phi <- function(dist_nu, sigma) {
    .Call(`_densityratio_make_Phi`, dist_nu, sigma)
}

make_phibar <- function(dist_de, sigma) {
    .Call(`_densityratio_make_phibar`, dist_de, sigma)
}

kliep_compute_alpha <- function(Phi, phibar, phibar_corr, epsilon, nepsilon, maxit, progressbar) {
    .Call(`_densityratio_kliep_compute_alpha`, Phi, phibar, phibar_corr, epsilon, nepsilon, maxit, progressbar)
}

compute_kliep <- function(dist_nu, dist_de, sigma, epsilon, maxit, cv_ind, progressbar) {
    .Call(`_densityratio_compute_kliep`, dist_nu, dist_de, sigma, epsilon, maxit, cv_ind, progressbar)
}

kmm_unconstrained_alpha <- function(Kdn, Kdd, Kd, nnu, nde) {
    .Call(`_densityratio_kmm_unconstrained_alpha`, Kdn, Kdd, Kd, nnu, nde)
}

kmm_constrained_alpha <- function(Kdn, Kdd, Kd, nnu, nde, settings) {
    .Call(`_densityratio_kmm_constrained_alpha`, Kdn, Kdd, Kd, nnu, nde, settings)
}

kmm_cv_loss <- function(Kdn, Kdd, Kd, Kn, nfolds, cv_ind_nu, cv_ind_de, constrained, settings) {
    .Call(`_densityratio_kmm_cv_loss`, Kdn, Kdd, Kd, Kn, nfolds, cv_ind_nu, cv_ind_de, constrained, settings)
}

compute_kmm <- function(nu, de, ce, Dd, sigma, cv_ind_nu, cv_ind_de, parallel, nthreads, progressbar, constrained, settings) {
    .Call(`_densityratio_compute_kmm`, nu, de, ce, Dd, sigma, cv_ind_nu, cv_ind_de, parallel, nthreads, progressbar, constrained, settings)
}

make_UV <- function(U) {
    .Call(`_densityratio_make_UV`, U)
}

get_sigma_lhss <- function(dist, sigma, quantiles) {
    .Call(`_densityratio_get_sigma_lhss`, dist, sigma, quantiles)
}

lhss_compute_alpha <- function(nu, de, ce, symmetric, m, intercept, sigma, quantiles, lambda, maxit, progressbar) {
    .Call(`_densityratio_lhss_compute_alpha`, nu, de, ce, symmetric, m, intercept, sigma, quantiles, lambda, maxit, progressbar)
}

compute_psihat <- function(K, Evecs, Evals, maxJ, ncol) {
    .Call(`_densityratio_compute_psihat`, K, Evecs, Evals, maxJ, ncol)
}

spectral_cv_loss <- function(Knu, Kde, m, maxM, nfolds, cv_ind_nu, cv_ind_de, nthreads, parallel) {
    .Call(`_densityratio_spectral_cv_loss`, Knu, Kde, m, maxM, nfolds, cv_ind_nu, cv_ind_de, nthreads, parallel)
}

spectral_dre <- function(dist_nu, dist_de, m, sigma, cv_ind_nu, cv_ind_de, parallel, nthreads, progressbar) {
    .Call(`_densityratio_spectral_dre`, dist_nu, dist_de, m, sigma, cv_ind_nu, cv_ind_de, parallel, nthreads, progressbar)
}

ulsif_compute_alpha <- function(Hhat, hhat, lambda) {
    .Call(`_densityratio_ulsif_compute_alpha`, Hhat, hhat, lambda)
}

set_threads <- function(nthreads) {
    .Call(`_densityratio_set_threads`, nthreads)
}

compute_ulsif_loocv <- function(Hhat, hhat, lambda, nnu, nde, nmin, ncol, Knu_nmin, Kde_nmin) {
    .Call(`_densityratio_compute_ulsif_loocv`, Hhat, hhat, lambda, nnu, nde, nmin, ncol, Knu_nmin, Kde_nmin)
}

compute_ulsif <- function(dist_nu, dist_de, sigma, lambda, parallel, nthreads, progressbar) {
    .Call(`_densityratio_compute_ulsif`, dist_nu, dist_de, sigma, lambda, parallel, nthreads, progressbar)
}

