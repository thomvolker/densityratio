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
#' @importFrom stats median
#'
#' @examples
#' x <- rnorm(100) |> matrix(25, 4)
#' densityratio:::kernel_gaussian(x)
#' y <- rnorm(200) |> matrix(50, 4)
#' densityratio:::kernel_gaussian(x, y)
#' sigma <- 5
#' densityratio:::kernel_gaussian(x, y, sigma)


kernel_gaussian <- function(x, y = NULL, sigma = NULL) {

  if(!is.null(y)) {
    distance <- distXY(x, y)
    if (is.null(sigma)) {
      sigma <- median(distance) |> sqrt()
    }
  } else {
    dist_vec <- distX(x)
    distance <- matrix(0, nrow(x), nrow(x))
    distance[lower.tri(distance)] <- dist_vec
    distance <- distance + t(distance)

    if (is.null(sigma)) {
      sigma <- median(distance[lower.tri(distance)]) |> sqrt()
    }
  }

  exp(-distance/(2*sigma^2))
}
