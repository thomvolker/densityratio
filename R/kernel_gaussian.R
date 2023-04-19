#' Compute Gaussian Kernel Gram matrix
#'
#' @param distance Numeric distance matrix computed with \code{distance}
#' @param symmetric Logical indicating whether the distance matrix is symmetric
#' which can speed up calculations.
#' @param sigma \code{NULL} or a scalar. If \code{NULL}, sigma is the
#' median Euclidean interpoint distance.
#' @return The Gaussian kernel gram matrix \eqn{\bf{K}}.
#' @importFrom stats median
#'
#' @examples
#' x <- rnorm(100) |> matrix(25, 4)
#' densityratio:::distance(x, x, TRUE) |> densityratio:::kernel_gaussian(symmetric = TRUE)
#' y <- rnorm(200) |> matrix(50, 4)
#' densityratio:::distance(x, y) |> densityratio:::kernel_gaussian()
#' sigma <- 5
#' densityratio:::distance(x, y) |> densityratio:::kernel_gaussian(sigma)

kernel_gaussian <- function(distance, sigma = NULL, symmetric = FALSE) {

  if (!is.null(sigma)) {
    stopifnot(is.numeric(sigma), length(sigma) == 1)
  }
  if (is.null(sigma)) {
    sigma <- median_distance(distance, symmetric = symmetric)
  }

  exp(-distance / (2 * sigma * sigma))
}

median_distance <- function(distance, symmetric = FALSE) {
  if (symmetric) {
    sqrt(median(distance[lower.tri(distance)]))
  } else {
    sqrt(median(distance))
  }
}
