#' Predict function for density object
#'
#' @keywords internal
predict.density <- function(object, newdata, lambda = 1e-9, log = FALSE) {
  if (missing(newdata)) {
    if (isTRUE(log)) return(log(object$y)) else return(object$y)
  }
  if (isTRUE(log)) {
    res <- approx(object$x, log(object$y), newdata)$y
    res[is.na(res)] <- log(lambda)
  } else {
    res <- approx(object$x, object$y, newdata)$y
    res[is.na(res)] <- lambda
  }
  return(res)
}

#' Naive approach to density ratio estimation
#'
#' The naive approach creates separate kernel density estimates for
#' the numerator and the denominator samples, and then evaluates their
#' ratio for the denominator samples. For multivariate data, this
#' approach assumes the variables are independent (naive Bayes assumption).
#'
#' @param nu Numeric matrix with numerator samples
#' @param de Numeric matrix with denominator samples (must have the same
#' variables as `nu`)
#' @param n the number of equally spaced points at which the density is to be
#' estimated. When n > 512, it is rounded up to a power of 2 during the
#' calculations (as fft is used) and the final result is interpolated by
#' [stats::approx]. So it almost always makes sense to specify n as a power of
#' two.
#' @param ... further arguments passed to [stats::density]
#'
#' @return `naive` returns `rhat_de`, the estimated density ratio for
#' the denominator samples.
#'
#' @importFrom stats density approx predict
#'
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(200, 1, 2)
#'
#' naive(x, y)
#' naive(x, y, bw = 2)
#'
#' @export
naive <- function(nu, de, n = 2L^11, ...) {
  nu <- as.matrix(nu)
  de <- as.matrix(de)
  N <- nrow(nu)
  P <- ncol(nu)

  # work on log scale
  # log-densities
  ld_nu_de <- numeric(N)
  ld_de_de <- numeric(N)

  # naive-bayes assumption of independence:
  # just add log-densities together
  for (p in 1:P) {
    d_nu_p <- density(nu[,p], n = n, ...)
    d_de_p <- density(de[,p], n = n, ...)
    ld_nu_de <- ld_nu_de + predict(d_nu_p, de[,p], log = TRUE)
    ld_de_de <- ld_de_de + predict(d_de_p, de[,p], log = TRUE)
  }

  return(exp(ld_nu_de - ld_de_de))
}

#' Naive subspace density ratio estimation
#'
#' The naive subspace estimator first creates an m-dimensional representation
#' of the data using singular value decomposition, and then runs the [naive]
#' density ratio estimation procedure on the data projected on this subspace.
#' The SVD is computed using the denominator samples.
#'
#' @param nu Numeric matrix with numerator samples
#' @param de Numeric matrix with denominator samples (must have the same
#' variables as `nu`)
#' @param m The size (in number of features) of the subspace
#' @param n the number of equally spaced points at which the density is to be
#' estimated. When n > 512, it is rounded up to a power of 2 during the
#' calculations (as fft is used) and the final result is interpolated by
#' [stats::approx]. So it almost always makes sense to specify n as a power of
#' two.
#' @param ... further arguments passed to [stats::density]
#'
#' @examples
#' set.seed(456)
#' # create data that differs only on the first variable
#' N <- 100
#' P <- 3 # P-1 noise variables
#' X <- matrix(rnorm(N*P), N)
#' Y <- cbind(rnorm(N, 1, 2), matrix(rnorm(N*P-N), N))
#' Y <- Y[order(Y[,1]),] # order so we can make nice line plot
#'
#' # plot: true, naive, naive_subspace
#' plot(Y[,1], dnorm(Y[,1]) / dnorm(Y[,1], 1, 2), type = "l", ylab = "Density ratio")
#' lines(Y[,1], naive(X, Y), col = "lightblue")
#' lines(Y[,1], naive_subspace(X, Y, 1), col = "darkorange")
#'
#' @export
naive_subspace <- function(nu, de, m, n = 2L^11, ...) {
  nu <- as.matrix(nu)
  de <- as.matrix(de)
  N <- nrow(nu)
  P <- ncol(nu)

  if (m > P) stop("Subspace size must be smaller than number of variables!")

  # first, use svd to compute m-dimensional subspace
  # then, run naive()
  # project de to m-dimensional space
  de_centered <- scale(de, scale = FALSE)
  V <- svd(de_centered, nu = m, nv = m)$v
  de_proj <- de_centered %*% V

  nu_centered <- scale(nu, center = attr(de_centered, "scaled:center"), scale = FALSE)
  nu_proj <- nu_centered %*% V

  return(naive(nu_proj, de_proj, ...))
}
