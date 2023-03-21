#' Kullback-Leibler importance estimation procedure
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
#' @export
#'
#' @return \code{kliep} returns \code{rhat_de}, the estimated density ratio for
#' the denominator samples.
#'
#' @examples
#' x <- rnorm(100) |> matrix(100)
#' y <- rnorm(200, 1, 2) |> matrix(200)
#' kliep(x, y)
#' kliep(x, y, sigma = 2, lambda = 2)
#'

kliep <- function(nu, de, sigma = NULL, maxcenters = nrow(nu), centers = NULL,
                  eps = 0.01, maxiter = 100) {

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
  if (!is.numeric(maxcenters) | length(maxcenters) > 1) {
    stop("maxcenters must be a scalar value indicating how many centers are maximally accepted (for computational reasons)")
  }
  if (is.null(centers)) {
    if (maxcenters < nrow(nu)) {
      centers <- nu[sample(nrow(nu), maxcenters), ]
    } else {
      centers <- nu
    }
  } else {
    centers <- as.matrix(centers)
    if (!is.numeric(centers) | !ncol(centers)==ncol(nu)) {
      stop("If centers are provided, they must have the same variables as the numerator samples")
    }
  }
    Phi <- kernel_gaussian(nu, centers, sigma) |> t()
    phibar <- kernel_gaussian(de, centers, sigma) |> colMeans() |> matrix()
    t_old <- runif(ncol(Phi))
    s_old <- mean(log(Phi %*% t_old))

    conv <- FALSE
    iter <- 0
    while(!conv) {
      iter <- iter+1
      cat(paste0("\r Iteration:", iter))
      t_new <- t_old + eps*t(Phi) %*% (1 / Phi %*% t_old)
      t_new <- t_new + c(1 - t(phibar)%*%t_new)*phibar/c(t(phibar)%*%phibar)
      t_new <- pmax(0, t_new)
      t_new <- t_new / c(t(phibar)%*%t_new)
      s_new <- mean(log(Phi %*% t_new))
      if (s_new <= s_old | iter == maxiter) conv <- TRUE
      #if (sum(abs(t_new - t_old)) < 1e-7 | iter == 100) conv <- TRUE
      t_old <- t_new
      s_old <- s_new
    }
    t_new
}

# N <- 500
# x <- rnorm(N) |> matrix(N)
# y <- rnorm(N, 1, 2) |> matrix(N)
#
# plot(y, densityratio:::kernel_gaussian(y, x) %*% kliep(x, y, 1.4, 5),
#      xlim = c(-5, 5), ylim = c(0, 2))
#
# mod_out <- densratio::KLIEP(x, y, kernel_num = N)
# plot(x, mod_out$compute_density_ratio(y))
# Psi <- kernel_gaussian(x) |> t()
# psib <- rowMeans(kernel_gaussian(x, y))
#
# theta <- runif(ncol(Psi))
# eps <- 0.01
# one <- rep(1, 50)
#

# .linear_kliep(Phi, phibar, eps, maxiter) {
#   theta_old <- rep(1, ncol(Phi))
#   conv <- FALSE
#
# }
