#' Least-squares heterodistributional subspace search
#'
#' @param df_numerator \code{data.frame} with exclusively numeric variables with
#' the numerator samples
#' @param df_denominator \code{data.frame} with exclusively numeric variables
#' with the denominator samples (must have the same variables as
#' \code{df_denominator})
#' @param m Scalar indicating the dimensionality of the reduced subspace
#' @param intercept \code{logical} Indicating whether to include an intercept
#' term in the model. Defaults to \code{TRUE}.
#' @param scale \code{"numerator"}, \code{"denominator"}, or \code{NULL},
#' indicating whether to standardize each numeric variable according to the
#' numerator means and standard deviations, the denominator means and standard
#' deviations, or apply no standardization at all.
#' @param nsigma Integer indicating the number of sigma values (bandwidth
#' parameter of the Gaussian kernel gram matrix) to use in cross-validation.
#' @param sigma_quantile \code{NULL} or numeric vector with probabilities to
#' calculate the quantiles of the distance matrix to obtain sigma values. If
#' \code{NULL}, \code{nsigma} values between \code{0.05} and \code{0.95} are
#' used.
#' @param sigma \code{NULL} or a scalar value to determine the bandwidth of the
#' Gaussian kernel gram matrix. If \code{NULL}, \code{nsigma} values between
#' \code{0.05} and \code{0.95} are used.
#' @param nlambda Integer indicating the number of \code{lambda} values
#' (regularization parameter), by default, \code{lambda} is set to
#' \code{10^seq(3, -3, length.out = nlambda)}.
#' @param lambda \code{NULL} or numeric vector indicating the lambda values to
#' use in cross-validation
#' @param ncenters Maximum number of Gaussian centers in the kernel gram
#' matrix. Defaults to all numerator samples.
#' @param centers Numeric matrix with the same variables as \code{nu} and
#' \code{de} that are used as Gaussian centers in the kernel Gram matrix. By
#' default, the matrix \code{nu} is used as the matrix with Gaussian centers.
#' @param maxit Maximum number of iterations in the updating scheme.
#' @param progressbar Logical indicating whether or not to display a progressbar.
#' @export
#' @return \code{lhss}-object, containing all information to calculate the
#' density ratio using optimal sigma, optimal lambda and optimal weights.
#'
#' @return \code{lhss} returns \code{rhat}, the estimated density ratio.
#' @example inst/examples/naive-example.R



lhss <- function(df_numerator, df_denominator, m = NULL, intercept = TRUE,
                 scale = "numerator", nsigma = 10, sigma_quantile = NULL,
                 sigma = NULL, nlambda = 10, lambda = NULL, ncenters = 200,
                 centers = NULL, maxit = 200, progressbar = TRUE) {

  cl <- match.call()
  nu <- check.datatype(df_numerator)
  de <- check.datatype(df_denominator)

  check.variables(nu, de, centers)

  df_centers <- check.centers(rbind(nu, de), centers, ncenters)
  dat <- check.dataform(nu, de, df_centers, is.null(centers), NULL, scale)

  p <- ncol(dat$nu)

  symmetric      <- check.symmetric(dat$nu, dat$ce)
  sigma_quantile <- check.sigma_quantile.lhss(nsigma, sigma, sigma_quantile)
  sigma          <- if (is.null(sigma_quantile)) sigma else sigma_quantile
  is_quantile    <- !is.null(sigma_quantile)
  lambda         <- check.lambda(nlambda, lambda)
  m              <- check.subspace(m, p)
  maxit          <- check.maxit(maxit)

  res <- lhss_compute_alpha(dat$nu, dat$de, dat$ce, symmetric, m, intercept, sigma,
                            is_quantile, lambda, maxit, progressbar)

  colnames(res$loocv) <- names(lambda) <- paste0("lambda", 1:length(lambda))
  rownames(res$loocv) <- names(sigma) <- paste0("sigma", 1:length(sigma))
  dimnames(res$sigmaopt) <- dimnames(res$loocv)
  dimnames(res$alpha) <- list(NULL, names(sigma), names(lambda))
  U <- array(res$Uopt, dim = c(p, m, length(sigma), length(lambda)),
             dimnames = list(NULL, NULL, names(sigma), names(lambda)))
  min_score <- arrayInd(which.min(res$loocv), dim(res$loocv))

  out <- list(
    df_numerator = df_numerator,
    df_denominator = df_denominator,
    alpha = res$alpha,
    cv_score = res$loocv,
    scale = scale,
    sigma = res$sigmaopt,
    sigma_quantiles = sigma_quantile,
    lambda = lambda,
    U = U,
    m = m,
    centers = df_centers,
    model_matrices = dat,
    alpha_opt = res$alpha[, min_score[1], min_score[2]],
    lambda_opt = lambda[min_score[2]],
    sigma_opt = res$sigmaopt[min_score[1], min_score[2]],
    U_opt = U[,,min_score[1], min_score[2], drop = FALSE],
    call = cl
  )
  class(out) <- "lhss"
  out
}
