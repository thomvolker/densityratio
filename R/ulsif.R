#' Unconstrained least-squares importance fitting
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
#' @examples
#' set.seed(1)
#' x <- rnorm(100) |> matrix(100)
#' y <- rnorm(200, 1, 2) |> matrix(200)
#' ulsif(x, y)
#' ulsif(x, y, sigma = 2, lambda = 2)

ulsif <- function(df_numerator, df_denominator, nsigma = 10, sigma_quantile = NULL, sigma = NULL,
                  nlambda = 20, lambda = NULL, ncenters = 200,
                  centers = NULL, parallel = FALSE, nthreads = NULL,
                  progressbar = TRUE) {

  cl <- match.call()
  nu <- as.matrix(df_numerator)
  de <- as.matrix(df_denominator)

  check.dataform(nu, de)
  centers   <- check.centers(nu, centers, ncenters)
  symmetric <- check.symmetric(nu, centers)
  parallel  <- check.parallel(parallel, nthreads, sigma, lambda)
  nthreads  <- check.threads(parallel, nthreads)

  dist_nu <- distance(nu, centers, symmetric)
  dist_de <- distance(de, centers)

  sigma  <- check.sigma(nsigma, sigma_quantile, sigma, dist_nu)
  lambda <- check.lambda(nlambda, lambda)

  res <- compute_ulsif(dist_nu, dist_de, sigma, lambda, parallel, nthreads, progressbar)
  min_score <- which.min(res$loocv_score) - 1
  alpha_min <- c(lambda =  min_score %% length(lambda) + 1,
                 sigma = min_score %/% length(lambda) + 1)

  out <- list(
    df_numerator = df_numerator,
    df_denominator = df_denominator,
    alpha = res$alpha,
    cv_score = res$loocv_score,
    sigma = sigma,
    lambda = lambda,
    centers = centers,
    alpha_opt = res$alpha[, alpha_min["lambda"], alpha_min["sigma"]],
    lambda_opt = lambda[alpha_min["lambda"]],
    sigma_opt = sigma[alpha_min["sigma"]],
    call = cl
  )
  class(out) <- "ulsif"
  out
}
