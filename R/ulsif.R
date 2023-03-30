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
#' @param maxcenters Maximum number of Gaussian centers in the kernel gram
#' matrix. Defaults to all numerator samples.
#' @param centers Numeric matrix with the same variables as \code{nu} and
#' \code{de} that are used as Gaussian centers in the kernel Gram matrix. By
#' default, the matrix \code{nu} is used as the matrix with Gaussian centers.
#' @export
#'
#' @return \code{ulsif} returns \code{rhat}, the estimated density ratio.
#'
#' @examples
#' x <- rnorm(100) |> matrix(100)
#' y <- rnorm(200, 1, 2) |> matrix(200)
#' ulsif(x, y)
#' ulsif(x, y, sigma = 2, lambda = 2)
#'

ulsif <- function(nu, de, sigma = NULL, lambda = NULL, maxcenters = nrow(nu),
                  centers = NULL) {

  nu <- as.matrix(nu)
  de <- as.matrix(de)

  n_nu <- nrow(nu)
  n_de <- nrow(de)
  p    <- ncol(nu)

  if (!is.numeric(de) | !is.numeric(nu) | !p == ncol(de)) {
    stop("Arguments de and nu must be matrices with the same variables")
  }
  if (!is.null(sigma)) {
    if (!is.numeric(sigma) | length(sigma) > 1) {
      stop("If provided, sigma must be a scalar")
    }
  }
  if (!is.null(lambda)) {
    if (!is.numeric(lambda) | length(lambda) > 1) {
      stop("Only scalar values for lambda are currently accepted")
    }
  } else {
    lambda <- sqrt(n_nu + n_de)
  }

  if (!is.numeric(maxcenters) | length(maxcenters) > 1) {
    stop("maxcenters must be a scalar value indicating how many centers are maximally accepted (for computational reasons)")
  }
  if (is.null(centers)) {
    if (maxcenters < nrow(nu)) {
      centers <- nu[sample(n_nu, maxcenters), ]
    } else {
      centers <- nu
    }
  } else {
    centers <- as.matrix(centers)
    if (!is.numeric(centers) | ! p == ncol(centers)) {
      stop("If centers are provided, they must have the same variables as the numerator samples")
    }
  }

 phi_nu <- kernel_gaussian(nu, centers, sigma)
 phi_de <- kernel_gaussian(de, centers, sigma)
 Hhat   <- crossprod(phi_de) / n_de
 one    <- diag(1, ncol(Hhat))
 hhat   <- colMeans(phi_nu)

 theta <- solve(Hhat + lambda * one) %*% hhat
 theta
}
