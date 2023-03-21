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
#' @importFrom quadprog solve.QP
#' @export
#'
#' @return \code{kmm} returns \code{rhat_de}, the estimated density ratio for
#' the denominator samples.
#'
#' @examples
#' x <- rnorm(100) |> matrix(100)
#' y <- rnorm(200, 1, 2) |> matrix(200)
#' kmm(x, y)
#' kmm(x, y, sigma = 2, lambda = 2)
#'

ulsif <- function(nu, de, sigma = NULL, lambda = NULL, maxcenters = nrow(nu)) {

  nu <- as.matrix(nu)
  de <- as.matrix(de)

  if (!is.numeric(de) | !is.numeric(nu) | !ncol(de) == ncol(nu)) {
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
    lambda <- sqrt(nrow(nu) + nrow(de))
  }
#  if ()

#  Kdede <- kernel_gaussian()
}
