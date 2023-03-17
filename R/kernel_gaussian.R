#' Compute Gaussian Kernel Gram matrix
#'
#' @param x Matrix to compute Gaussian kernel gram matrix with
#' @param y \code{NULL} or a matrix with the same variables as \code{x}.
#' If \code{NULL}, the kernel gram matrix of \code{x} with itself is
#' computed, otherwise, the kernel gram matrix of \code{x} with
#' \code{y} is computed.
#' @param sigma \code{NULL} or a scalar. If \code{NULL}, sigma is the
#' median Euclidean interpoint distance.
#' @return The Gaussian kernel gram matrix \eqn{\bf{K}}.
#' @export kernel_gaussian
#'
#' @examples
#' x <- rnorm(100) |> matrix(25, 4)
#' kernel_gaussian(x)
#' y <- rnorm(200) |> matrix(50, 4)
#' kernel_gaussian(x, y)
#' sigma <- 5
#' kernel_gaussian(x, y, sigma)


kernel_gaussian <- function(x, y = NULL, sigma = NULL) {

  x <- as.matrix(x)

  if (!is.numeric(x)[1]) stop("x must be a numeric matrix")

  if(!is.null(y)) {
    y <- as.matrix(y)
    if ((!is.numeric(y)[1]) | (!ncol(x) == ncol(y))) {
      stop("y must be a numeric matrix with the same number of columns as x")
    }
    distance <- distXY(x, y, nrow(x), nrow(y), ncol(x))
  } else {
    distance <- matrix(0, nrow(x), nrow(x))
    distance[lower.tri(distance)] <- distX(x, nrow(x), ncol(x))
    distance <- distance + t(distance)
  }
  if (is.null(sigma)) {
    sigma <- stats::median(distance) |> sqrt()
  } else if (!is.numeric(sigma) | length(sigma) > 1) {
    stop("sigma, if supplied, must be a scalar")
  }
  exp(-distance/(2*sigma^2))
}
