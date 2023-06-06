#' Unconstrained least-squares importance fitting
#'
#' @param nu Numeric matrix with numerator samples
#' @param de Numeric matrix with denominator samples (must have the same
#' variables as \code{nu})
#' @param sigma \code{NULL} or a scalar value to determine the bandwidth of the
#' Gaussian kernel gram matrix. If \code{NULL}, sigma is the median Euclidean
#' interpoint distance.
#' @param lambda \code{NULL} or a scalar value to determine the regularization
#' imposed on the Gaussian kernel gram matrix of the denominator samples. If
#' \code{NULL}, \code{lambda} is chosen to be \eqn{\sqrt{N}}.
#' @param ncenters Maximum number of Gaussian centers in the kernel gram
#' matrix. Defaults to all numerator samples.
#' @param centers Numeric matrix with the same variables as \code{nu} and
#' \code{de} that are used as Gaussian centers in the kernel Gram matrix. By
#' default, the matrix \code{nu} is used as the matrix with Gaussian centers.
#' @param parallel Logical argument indicating whether the density ratio
#' parameters should be calculated in parallel for different lambda values.
#' @param nthreads Scalar value indicating the number of threads to use for
#' parallel processing.
#' @export
#'
#' @return \code{ulsif} returns \code{rhat}, the estimated density ratio.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(100) |> matrix(100)
#' y <- rnorm(200, 1, 2) |> matrix(200)
#' ulsif(x, y)
#' ulsif(x, y, sigma = 2, lambda = 2)

ulsif <- function(nu, de, nsigma = 10, sigma_quantile = NULL, sigma = NULL,
                  nlambda = 20, lambda = NULL, ncenters = 200,
                  centers = NULL, parallel = FALSE, nthreads = NULL,
                  progressbar = TRUE) {

  nu <- as.matrix(nu)
  de <- as.matrix(de)

  n_nu <- nrow(nu)
  n_de <- nrow(de)
  p    <- ncol(nu)

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
    alpha = res$alpha,
    loocv_score = res$loocv_score,
    sigma = sigma,
    lambda = lambda,
    centers = centers,
    alpha_min = res$alpha[, alpha_min["lambda"], alpha_min["sigma"]],
    lambda_min = lambda[alpha_min["lambda"]],
    sigma_min = sigma[alpha_min["sigma"]]
  )
  class(out) <- "dratio"
  out
}
