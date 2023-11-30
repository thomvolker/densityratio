#' Naive density ratio estimation
#'
#' The naive approach creates separate kernel density estimates for
#' the numerator and the denominator samples, and then evaluates their
#' ratio for the denominator samples. For multivariate data, this
#' approach assumes the variables are independent (naive Bayes assumption).
#'
#' @param df_numerator \code{data.frame} with exclusively numeric variables with
#' the numerator samples
#' @param df_denominator \code{data.frame} with exclusively numeric variables
#' with the denominator samples (must have the same variables as
#' \code{df_denominator})
#' @param n \code{integer} the number of equally spaced points at which the density is
#' estimated. When n > 512, it is rounded up to a power of 2 during the
#' calculations (as fft is used) and the final result is interpolated by
#' [stats::approx]. So it makes sense to specify n as a power of two.
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
naive <- function(df_numerator, df_denominator, n = 2L^11, ...) {
  cl <- match.call()
  nu <- as.matrix(df_numerator)
  de <- as.matrix(df_denominator)
  P <- ncol(nu)
  check.dataform(nu, de)

  # naive-bayes like assumption of independence:
  # compute and store densities for each column
  d_nu <- lapply(1:P, \(p) density(nu[,p], n = n, ...))
  d_de <- lapply(1:P, \(p) density(de[,p], n = n, ...))

  # return object
  out <- list(
    df_numerator = df_numerator,
    df_denominator = df_denominator,
    density_numerator = d_nu,
    density_denominator = d_de,
    call = cl
  )
  class(out) <- c("naivedensityratio")
  return(out)
}

