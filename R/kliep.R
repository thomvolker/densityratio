#' Kullback-Leibler importance estimation procedure
#'
#' @param df_numerator \code{data.frame} with exclusively numeric variables with
#' the numerator samples
#' @param df_denominator \code{data.frame} with exclusively numeric variables
#' with the denominator samples (must have the same variables as
#' \code{df_denominator})
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
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(100) |> matrix(100)
#' y <- rnorm(200, 1, 2) |> matrix(200)
#' kliep(x, y)
#' kliep(x, y, nsigma = 20, ncenters = 100, nfold = 20, epsilon = 10^{3:-5}, maxit = 1000)
#'

kliep <- function(df_numerator, df_denominator, nsigma = 10, sigma_quantile = NULL, sigma = NULL,
                  ncenters = 200, centers = NULL, cv = TRUE, nfold = NULL,
                  epsilon = NULL, maxit = 2000, progressbar = TRUE) {

  cl <- match.call()
  nu <- as.matrix(df_numerator)
  de <- as.matrix(df_denominator)

  nnu <- nrow(nu)

  check.dataform(nu, de)
  centers   <- check.centers(nu, centers, ncenters)
  symmetric <- check.symmetric(nu, centers)

  dist_nu <- distance(nu, centers, symmetric)
  dist_de <- distance(de, centers)

  sigma   <- check.sigma(nsigma, sigma_quantile, sigma, dist_nu)
  epsilon <- check.epsilon(epsilon)
  maxit   <- check.maxit(maxit)
  cv_ind  <- check.nfold(cv, nfold, sigma, nnu)

  res <- compute_kliep(dist_nu, dist_de, sigma, epsilon, maxit, cv_ind, progressbar)
  rownames(res$cv_score) <- names(sigma) <- paste0("sigma", 1:length(sigma))
  dimnames(res$alpha) <- list(NULL, names(sigma))


  out <- list(
    df_numerator = df_numerator,
    df_denominator = df_denominator,
    alpha = res$alpha,
    cv_score = switch(cv, res$cv_score, NULL),
    sigma = sigma,
    centers = centers,
    nfold = switch(cv, max(cv_ind) + 1, NULL),
    epsilon = epsilon,
    maxit = maxit,
    alpha_opt = switch(cv, res$alpha[, which.max(res$cv_score), drop = FALSE], NULL),
    sigma_opt = switch(cv, sigma[which.max(res$cv_score)], NULL),
    call = cl
  )
  class(out) <- "kliep"
  out
}
