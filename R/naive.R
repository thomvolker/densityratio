#' Naive density ratio estimation
#'
#' The naive approach creates separate kernel density estimates for
#' the numerator and the denominator samples, and then evaluates their
#' ratio for the denominator samples. For multivariate data, the density ratio
#' is computed after a orthogonal linear transformation, such that the new
#' variables can be treated as independent. To reduce the dimensionality of
#' the PCA solution, one can set the number of components by setting the
#' \code{m} parameter to an integer value smaller than the number of variables.
#'
#' @param df_numerator \code{data.frame} with exclusively numeric variables with
#' the numerator samples
#' @param df_denominator \code{data.frame} with exclusively numeric variables
#' with the denominator samples (must have the same variables as
#' \code{df_denominator})
#' @param bw the smoothing bandwidth to be used. See [stats::density] for more
#' information.
#' @param kernel the kernel to be used. See [stats::density] for more
#' information.
#' @param m \code{integer} Optional parameter to reduce the dimensionality of
#' the data in multivariate density ratio estimation problems. If missing,
#' the number of variables in the data is used. If set to an integer value
#' smaller than the number of variables, the first \code{m} principal
#' components are used to estimate the density ratio. If set to \code{NULL},
#' the square root of the number of variables is used (for consistency with
#' other methods).
#' @param n \code{integer} the number of equally spaced points at which the
#' density is to be estimated. When n > 512, it is rounded up to a power of 2
#' during the calculations (as fast Fourier transform is used) and the final
#' result is interpolated by [stats::approx]. So it makes sense to specify n as
#' a power' of two.
#' @param ... further arguments passed to [stats::density]
#'
#' @return `naivedensityratio` object
#'
#' @importFrom stats density
#'
#' @seealso [stats::density()]
#'
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(200, 1, 2)
#'
#' naive(x, y)
#' naive(x, y, bw = 2)
#'
#' @export
naive <- function(df_numerator, df_denominator, m, bw = "SJ",
                  kernel = "gaussian", n = 2L^11, ...) {
  cl <- match.call()

  nu <- check.datatype(df_numerator)
  de <- check.datatype(df_denominator)

  check.variables(nu, de)
  dat <- check.dataform(nu, de, NULL, TRUE, NULL, NULL)

  P <- ncol(dat$nu)
  M <- ifelse(missing(m), P, m)
  M <- check.subspace(M, P)

  # PCA to remove linear relations between variables
  pca <- stats::prcomp(dat$nu, center = TRUE, scale. = TRUE, rank = M)
  nu_proj <- pca$x |> asplit(2)
  de_proj <- stats::predict(pca, newdata = dat$de) |> asplit(2)

  # now assume independence and calculate the density for each component
  # separately
  d_nu <- lapply(nu_proj, \(x) density(x, bw = bw, kernel = kernel, n = n, ...))
  d_de <- lapply(de_proj, \(x) density(x, bw = bw, kernel = kernel, n = n, ...))

  # return object
  out <- list(
    df_numerator = df_numerator,
    df_denominator = df_denominator,
    fit = pca,
    model_matrices = list(nu = nu, de = de),
    density_numerator = d_nu,
    density_denominator = d_de,
    call = cl
  )
  class(out) <- c("naivedensityratio")
  return(out)
}

