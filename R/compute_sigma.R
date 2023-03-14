
#' Compute median interpoint distance sigma
#'
#' @param X A data matrix containing the observations to compute the median interpoint distance on
#' @param nsamples The maximum number of samples that is considered
#'
#' @return A scalar \eqn{\sigma} (the median interpoint distance)
#' @export compute_sigma
#'
#' @examples
#' x <- rnorm(100) |> matrix(25, 4)
#' compute_sigma(x)


compute_sigma <- function(X, nsamples = nrow(X)) {
  d <- nrow(X)
  if (d > nsamples) { # weird step, it seems more logical to sample from X
    Xmed <- X[1:nsamples, ]
    d <- nsamples
  } else {
    Xmed <- X
  }
  G <- rowSums(Xmed*Xmed)
  Q <- rep(G, d) |> matrix(d, d)

  dist <- Q + t(Q) - 2*Xmed%*%t(Xmed)
  sigma <- stats::median(dist[upper.tri(dist)] / 2) |> sqrt()

  sigma
}
