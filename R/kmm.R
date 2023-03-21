#' Kernel mean matching approach to density ratio estimation
#'
#' @param nu Numeric matrix with numerator samples
#' @param de Numeric matrix with denominator samples (must have the same
#' variables as \code{nu})
#' @param method Character string containing the method used for kernel mean
#' matching. Currently, \code{method = "unconstrained"} and
#' \code{method = "constrained"} are supported.
#' @param sigma \code{NULL} or a scalar value to determine the bandwidth of the
#' Gaussian kernel gram matrix. If \code{NULL}, \code{sigma} is the median Euclidean
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

kmm <- function(nu, de, method = "unconstrained", sigma = NULL, lambda = NULL) {

  nu <- as.matrix(nu)
  de <- as.matrix(de)

  methods <- c("unconstrained", "constrained")

  if (length(method) > 1 | !method %in% methods) {
    stop(paste("Method must be a single method: constrained or unconstrained"))
  }
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

  Kdede <- kernel_gaussian(de, sigma = sigma)
  Kdenu <- kernel_gaussian(de, nu, sigma = sigma)

  if (method == "unconstrained") {
    rhat_de <- .kmm_unconstrained(Kdede, Kdenu, lambda)
  }
  if (method == "constrained") {
    rhat_de <- .kmm_constrained(Kdede, Kdenu, lambda)
  }
  list(rhat_de = rhat_de)

}

.kmm_unconstrained <- function(Kdede, Kdenu, lambda) {
  nnu <- ncol(Kdenu)
  nde <- nrow(Kdede)
  one <- rep(1, nnu)
  nde/nnu * solve(Kdede + lambda*diag(nde)) %*% Kdenu %*% one
}

.kmm_constrained <- function(Kdede, Kdenu, lambda) {

  # We impose that the density ratio multiplied with the denominator density
  # should integrate to 1 (up to some error eps), that the density ratio is
  # positive everywhere, and that the density ratio is constrained to be
  # smaller than some large positive scalar (here chosen to be 10000)

  nnu <- ncol(Kdenu)
  nde <- nrow(Kdede)
  one <- rep(1, nnu)
  # Constraints that make sure that the density ratio with denominator density
  # integrates to one, and the second parts sets the boundaries
  A <- cbind(rep(-1, nde),
             rep(1,  nde),
             diag(1,  nde, nde),
             diag(-1, nde, nde))

  # Bounds with respect to the integration
  eps <- (sqrt(nde) - 1) / (sqrt(nde))

  # Second part of constraints, t(A) %*% r >= b
  b <- c(-nde * (eps + 1),
         nde * (eps - 1),
         rep(0, nde),
         -rep(1000, nde))

  # Solve the optimization problem
  solve.QP(Dmat = Kdede + lambda*diag(nde),
           dvec = Kdenu %*% one,
           Amat = A,
           bvec = b)$solution
}


