#' Kernel mean matching approach to density ratio estimation
#'
#' @param df_numerator \code{data.frame} with exclusively numeric variables with
#' the numerator samples
#' @param df_denominator \code{data.frame} with exclusively numeric variables
#' with the denominator samples (must have the same variables as
#' \code{df_denominator})
#' @param scale \code{"numerator"}, \code{"denominator"}, or \code{FALSE},
#' indicating whether to standardize each numeric variable according to the
#' numerator means and standard deviations, the denominator means and standard
#' deviations, or apply no standardization at all.
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

kmm <- function(df_numerator, df_denominator, scale = "numerator",
                method = "unconstrained", sigma = NULL, lambda = NULL) {

  cl <- match.call()

  nu <- check.datatype(df_numerator)
  de <- check.datatype(df_denominator)

  check.variables(nu, de)

  dat <- check.dataform(nu, de, de, TRUE, NULL, scale)

  methods <- c("unconstrained", "constrained")

  # TODO: checks not implemented because currently cross-validation is not an option
  # in Kernel mean matching. This is something for the future.

  if (length(method) > 1 | !method %in% methods) {
    stop(paste("Method must be a single method: constrained or unconstrained"))
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
    lambda <- sqrt(nrow(dat$nu) + nrow(dat$de))
  }

  distdede <- distance(dat$de, dat$de, FALSE)
  distdenu <- distance(dat$de, dat$nu, FALSE)
  if (is.null(sigma)) {
    sigma <- stats::median(distdede[distdede > 0])
  }


  Kdede <- kernel_gaussian(distdede, sigma = sigma)
  Kdenu <- kernel_gaussian(distdenu, sigma = sigma)

  if (method == "unconstrained") {
    rhat_de <- .kmm_unconstrained(Kdede, Kdenu, lambda)
  }
  if (method == "constrained") {
    rhat_de <- .kmm_constrained(Kdede, Kdenu, lambda)
  }

  out <- list(rhat_de = rhat_de,
              sigma = sigma,
              lambda = lambda,
              call = cl)
  class(out) <- "kmm"

  out

}

.kmm_unconstrained <- function(Kdede, Kdenu, lambda) {
  nnu <- ncol(Kdenu)
  nde <- nrow(Kdede)
  nde/nnu * solve(Kdede + lambda*diag(nde), Kdenu %*% rep(1, nnu))
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
  quadprog::solve.QP(Dmat = Kdede + lambda*diag(nde),
                     dvec = Kdenu %*% one,
                     Amat = A,
                     bvec = b)$solution
}


