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

ulsif <- function(nu, de, sigma = NULL, lambda = NULL, ncenters = nrow(nu),
                  centers = NULL, parallel = FALSE, nthreads = NULL) {

  nu <- as.matrix(nu)
  de <- as.matrix(de)

  n_nu <- nrow(nu)
  n_de <- nrow(de)
  p    <- ncol(nu)

  check.dataform(nu, de)
  # TODO: Expand default sigma options (i.e., allow different default options (median distance, (LOO)CV, etc.))
  check.sigma(sigma)
  check.lambda(lambda)
  centers   <- check.centers(nu, centers, ncenters)
  symmetric <- check.symmetric(centers, ncenters)
  parallel  <- check.parallel(parallel, nthreads, sigma, lambda)
  nthreads  <- check.threads(parallel, nthreads)

  phi_nu <- kernel_gaussian(distance(nu, centers, symmetric), sigma, symmetric)
  phi_de <- kernel_gaussian(distance(de, centers), sigma)
  Hhat   <- crossprod(phi_de) / n_de
  hhat   <- colMeans(phi_nu)

  if (is.null(lambda)) {
    lambda <- 1
    # TODO: Better default, lambda_max in glmnet ||XTY||_\infty (max value of XTY) see hastie, tibs, tibs 2020
  }

  theta <- compute_ulsif(Hhat, hhat, lambda, parallel, nthreads)
  theta
}
