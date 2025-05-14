#' Kullback-Leibler importance estimation procedure
#'
#' @param df_numerator \code{data.frame} with exclusively numeric variables with
#' the numerator samples
#' @param df_denominator \code{data.frame} with exclusively numeric variables
#' with the denominator samples (must have the same variables as
#' \code{df_denominator})
#' @param scale \code{"numerator"}, \code{"denominator"}, or \code{NULL},
#' indicating whether to standardize each numeric variable according to the
#' numerator means and standard deviations, the denominator means and standard
#' deviations, or apply no standardization at all.
#' @param nsigma Integer indicating the number of sigma values (bandwidth
#' parameter of the Gaussian kernel gram matrix) to use in cross-validation.
#' @param sigma_quantile \code{NULL} or numeric vector with probabilities to
#' calculate the quantiles of the distance matrix to obtain sigma values. If
#' \code{NULL}, \code{nsigma} values between \code{0.25} and \code{0.75} are
#' used.
#' @param sigma \code{NULL} or a scalar value to determine the bandwidth of the
#' Gaussian kernel gram matrix. If \code{NULL}, \code{nsigma} values between
#' \code{0.25} and \code{0.75} are used.
#' @param ncenters Maximum number of Gaussian centers in the kernel gram
#' matrix. Defaults to all numerator samples.
#' @param centers Option to specify the Gaussian samples manually.
#' @param cv Logical indicating whether or not to do cross-validation
#' @param nfold Number of cross-validation folds used in order to calculate the
#' optimal \code{sigma} value (default is 5-fold cv).
#' @param epsilon Numeric scalar or vector with the learning rate for the
#' gradient-ascent procedure. If a vector, all values are used as the learning
#' rate. By default, \code{10^{1:-5}} is used.
#' @param maxit Maximum number of iterations for the optimization scheme.
#' @param progressbar Logical indicating whether or not to display a progressbar.
#' @export
#'
#' @return \code{kliep}-object, containing all information to calculate the
#' density ratio using optimal sigma and optimal weights.
#' @references Sugiyama, M., Suzuki, T., Nakajima, S., Kashima, H., Von Bünau, P.,
#' & Kawanabe, M. (2008). Direct importance estimation for covariate shift
#' adaptation. <i>Annals of the Institute of Statistical Mathematics</i>
#' 60:699-746. Doi: https://doi.org/10.1007/s10463-008-0197-x.
#'
#' @example inst/examples/kliep-example.R

kliep <- function(df_numerator, df_denominator, scale = "numerator", nsigma = 10,
                  sigma_quantile = NULL, sigma = NULL, ncenters = 200,
                  centers = NULL, cv = TRUE, nfold = 5, epsilon = NULL,
                  maxit = 5000, progressbar = TRUE) {
  cl <- match.call()
  nu <- check.datatype(df_numerator)
  de <- check.datatype(df_denominator)

  check.variables(nu, de, centers)

  df_centers <- check.centers(rbind(nu, de), centers, ncenters)
  dat <- check.dataform(nu, de, df_centers, is.null(centers), NULL, scale)

  nnu <- nrow(dat$nu)

  dist_nu <- distance(dat$nu, dat$ce)
  dist_de <- distance(dat$de, dat$ce)

  sigma <- check.sigma(nsigma, sigma_quantile, sigma, dist_nu)
  epsilon <- check.epsilon(epsilon)
  maxit <- check.maxit(maxit)
  cv_ind <- check.nfold(cv, nfold, nnu)

  res <- compute_kliep(dist_nu, dist_de, sigma, epsilon, maxit, cv_ind, progressbar)
  rownames(res$cv_score) <- names(sigma) <- paste0("sigma", 1:length(sigma))
  dimnames(res$alpha) <- list(NULL, names(sigma))


  out <- list(
    df_numerator = df_numerator,
    df_denominator = df_denominator,
    alpha = res$alpha,
    cv_score = switch(cv,
      res$cv_score,
      NULL
    ),
    scale = scale,
    sigma = sigma,
    centers = df_centers,
    model_matrices = dat,
    nfold = switch(cv,
      max(cv_ind) + 1,
      NULL
    ),
    epsilon = epsilon,
    maxit = maxit,
    alpha_opt = switch(cv,
      res$alpha[, which.max(res$cv_score), drop = FALSE],
      NULL
    ),
    sigma_opt = switch(cv,
      sigma[which.max(res$cv_score)],
      NULL
    ),
    call = cl
  )
  class(out) <- "kliep"
  out
}
