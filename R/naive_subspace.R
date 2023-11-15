#' Naive subspace density ratio estimation
#'
#' The naive subspace estimator first creates an m-dimensional representation
#' of the data using singular value decomposition, and then runs the [naive]
#' density ratio estimation procedure on the data projected on this subspace.
#' The SVD is computed using the denominator samples.
#'
#' @param df_numerator \code{data.frame} with exclusively numeric variables with
#' the numerator samples
#' @param df_denominator \code{data.frame} with exclusively numeric variables
#' with the denominator samples (must have the same variables as
#' \code{df_denominator})
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
#' df_X <- data.frame(X)
#' df_Y <- data.frame(Y)
#'
#' # estimate
#' dr_naive <- naive(df_X, df_Y)
#' dr_subspace <- naivesubspace(df_X, df_Y, 1)
#'
#' # plot: true, naive, naive_subspace
#' df_new <- data.frame(
#'   X1 = seq(-4, 4, length.out = 100),
#'   X2 = seq(-4, 4, length.out = 100),
#'   X3 = seq(-4, 4, length.out = 100)
#' )
#' true_dr <- dnorm(df_new[,1]) / dnorm(df_new[,1], 1, 2)
#' plot(df_new[,1], true_dr, type = "l", ylab = "Density ratio", ylim = c(0, 3.2))
#' lines(df_new[,1], predict(dr_naive, df_new), col = "lightblue")
#' lines(df_new[,1], predict(dr_subspace, df_new), col = "darkorange")
#'
#' @export
naivesubspace <- function(df_numerator, df_denominator, m = NULL, n = 2L^11, ...) {
  cl <- match.call()
  nu <- as.matrix(df_numerator)
  de <- as.matrix(df_denominator)
  check.dataform(nu, de)

  m <- check.subspace(m, ncol(nu))

  # first, use svd to compute m-dimensional subspace
  de_centered <- scale(de, scale = FALSE)
  center <- attr(de_centered, "scaled:center")
  V <- svd(de_centered, nu = m, nv = m)$v
  de_proj <- de_centered %*% V

  nu_centered <- scale(nu, center = center, scale = FALSE)
  nu_proj <- nu_centered %*% V

  # then, perform naive density ratio estimation
  d_nu <- lapply(1:m, \(p) density(nu_proj[,p], n = n, ...))
  d_de <- lapply(1:m, \(p) density(de_proj[,p], n = n, ...))

  # return object
  out <- list(
    df_numerator = df_numerator,
    df_denominator = df_denominator,
    projection_matrix = V,
    subspace_dim = m,
    center = center,
    density_numerator = d_nu,
    density_denominator = d_de,
    call = cl
  )
  class(out) <- c("naivesubspacedensityratio")
  return(out)
}
