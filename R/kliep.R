#' Kullback-Leibler importance estimation procedure
#'
#' @param nu Numeric matrix with numerator samples
#' @param de Numeric matrix with denominator samples (must have the same
#' variables as \code{nu})
#' @param sigma \code{NULL} or a scalar value to determine the bandwidth of the
#' Gaussian kernel gram matrix. If \code{NULL}, sigma is the median Euclidean
#' interpoint distance.
#' @param ncenters Maximum number of Gaussian centers in the kernel gram
#' matrix. Defaults to all numerator samples.
#' @param centers Option to specify the Gaussian samples manually.
#' @param eps Learning rate for the optimization procedure (can be a vector)
#' @param maxit Maximum number of iterations for the optimization scheme.
#' @param printFlag Whether or not to print the progress of the optimization
#' scheme.
#' @export
#'
#' @return \code{kliep} returns \code{rhat_de}, the estimated density ratio for
#' the denominator samples.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(100) |> matrix(100)
#' y <- rnorm(200, 1, 2) |> matrix(200)
#' kliep(x, y)
#' kliep(x, y, sigma = 2)
#'

kliep <- function(nu, de, sigma = NULL, ncenters = nrow(nu), centers = NULL,
                  eps = 0.001, maxit = 100, printFlag = T) {

  nu <- as.matrix(nu)
  de <- as.matrix(de)

  n_nu <- nrow(nu)
  n_de <- nrow(de)
  p    <- ncol(nu)

  if (!is.numeric(de) | !is.numeric(nu) | ! p == ncol(de)) {
    stop("Arguments de and nu must be matrices with the same variables")
  }
  if (!is.null(sigma)) {
    if (!is.numeric(sigma) | length(sigma) > 1) {
      stop("If provided, sigma must be a scalar")
    }
  }
  if (!is.numeric(ncenters) | length(ncenters) > 1) {
    stop("ncenters must be a scalar value indicating how many centers are maximally accepted (for computational reasons)")
  }
  if (is.null(centers)) {
    if (ncenters < n_nu) {
      centers <- nu[sample(n_nu, ncenters), ]
    } else {
      centers <- nu
    }
  } else {
    centers <- as.matrix(centers)
    ncenters <- nrow(centers)
    if (!is.numeric(centers) | ! p == ncol(centers)) {
      stop("If centers are provided, they must have the same variables as the numerator samples")
    }
  }

  Phi <- kernel_gaussian(distance(nu, centers), sigma) |> t()
  phibar <- kernel_gaussian(distance(de, centers), sigma) |> colMeans()
  phibar_cp_phibar_inv <- phibar / c(crossprod(phibar))

  theta <- rep(1, ncenters)
  theta <- .impose_constraints(theta, phibar, phibar_cp_phibar_inv)
  score <- mean(log(crossprod(Phi, theta)))

  conv <- FALSE
  iter <- 0

  for (e in eps) {
    conv <- FALSE
    iter <- 0

    while (!conv) {
      iter <- iter+1
      if(printFlag) cat(paste0("\r Iteration: ", iter))
      t_temp  <- .gradient_ascent(theta, Phi, e)
      t_new   <- .impose_constraints(t_temp, phibar, phibar_cp_phibar_inv)
      s_new   <- mean(log(crossprod(Phi, t_new)))
      if (s_new - score <= 0 | iter == maxit) {
        conv <- TRUE
      } else {
        score <- s_new
        theta <- t_new
      }
    }
  }
  theta
}

.gradient_ascent <- function(theta, Phi, eps) {
  theta + eps * crossprod(Phi, 1 / Phi %*% theta)
}

.impose_constraints <- function(theta, phibar, phibar_cp_phibar_inv) {
  theta <- theta + phibar_cp_phibar_inv * c(1 - crossprod(phibar, theta))
  theta <- pmax(0, theta)
  theta %*% (1 / crossprod(phibar, theta))
}
