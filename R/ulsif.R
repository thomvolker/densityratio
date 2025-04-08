#' Unconstrained least-squares importance fitting
#'
#' @param df_numerator \code{data.frame} with exclusively numeric variables with
#' the numerator samples
#' @param df_denominator \code{data.frame} with exclusively numeric variables
#' with the denominator samples (must have the same variables as
#' \code{df_denominator})
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
#' @param centers \code{NULL} or numeric matrix with the same dimensions as the
#' data, indicating the centers for the Gaussian kernel gram matrix.
#' @param parallel logical indicating whether to use parallel processing in the
#' cross-validation scheme.
#' @param nthreads \code{NULL} or integer indicating the number of threads to
#' use for parallel processing. If parallel processing is enabled, it defaults
#' to the number of available threads minus one.
#' @param progressbar Logical indicating whether or not to display a progressbar.
#' @export
#'
#' @return \code{ulsif}-object, containing all information to calculate the
#' density ratio using optimal sigma and optimal weights.
#'
#' @example inst/examples/ulsif-example.R

ulsif <- function(df_numerator, df_denominator, intercept = TRUE, scale = "numerator",
                  nsigma = 10, sigma_quantile = NULL, sigma = NULL, nlambda = 20,
                  lambda = NULL, ncenters = 200, centers = NULL,
                  parallel = FALSE, nthreads = NULL, progressbar = TRUE) {
  cl <- match.call()
  nu <- check.datatype(df_numerator)
  de <- check.datatype(df_denominator)

  check.variables(nu, de, centers)

  df_centers <- check.centers(rbind(nu, de), centers, ncenters)
  dat <- check.dataform(nu, de, df_centers, is.null(centers), NULL, scale)

  intercept <- check.intercept(intercept)

  dist_nu <- distance(dat$nu, dat$ce, intercept)
  dist_de <- distance(dat$de, dat$ce, intercept)

  sigma <- check.sigma(nsigma, sigma_quantile, sigma, dist_nu)
  lambda <- check.lambda(nlambda, lambda)

  parallel <- check.parallel(parallel, nthreads, lambda)
  nthreads <- check.threads(parallel, nthreads)

  res <- compute_ulsif(dist_nu, dist_de, sigma, lambda, parallel, nthreads, progressbar)
  loocv_scores <- res$loocv_score
  colnames(loocv_scores) <- names(lambda) <- paste0("lambda", 1:length(lambda))
  rownames(loocv_scores) <- names(sigma) <- paste0("sigma", 1:length(sigma))
  dimnames(res$alpha) <- list(NULL, names(sigma), names(lambda))
  min_score <- arrayInd(which.min(loocv_scores), dim(loocv_scores))

  out <- list(
    df_numerator = df_numerator,
    df_denominator = df_denominator,
    alpha = res$alpha,
    cv_score = loocv_scores,
    scale = scale,
    sigma = sigma,
    lambda = lambda,
    centers = df_centers,
    model_matrices = dat,
    alpha_opt = res$alpha[, min_score[1], min_score[2]],
    lambda_opt = lambda[min_score[2]],
    sigma_opt = sigma[min_score[1]],
    call = cl
  )
  class(out) <- "ulsif"
  out
}
