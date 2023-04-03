#' Predict function for density ratio estimation
#'
#' @keywords internal
predict.density <- function(object, newdata, lambda = 1e-9) {
  if (missing(newdata)) {
    return(object$y)
  }
  res <- approx(object$x, object$y, newdata)$y
  res[is.na(res)] <- lambda
  res
}

#' Naive density ratio estimation
#'
#' @param nu Numeric matrix with numerator samples
#' @param de Numeric matrix with denominator samples (must have the same
#' variables as \code{nu})
#' @param ... arguments passed to [stats::density]
#'
#' @export
naive <- function(nu, de, ...) {
  d_nu <- density(nu, ...)
  d_de <- density(de, ...)
  d_nu_de <- predict(d_nu, de)
  d_de_de <- predict(d_de, de)
  return(d_nu_de / d_de_de)
}

#' Naive density ratio estimation
#'
#' @param nu Numeric matrix with numerator samples
#' @param de Numeric matrix with denominator samples (must have the same
#' variables as \code{nu})
#' @param m The size (in number of features) of the subspace
#' @param ... arguments passed to [stats::density]
#'
#' @export
naive_subspace <- function(nu, de, m, ...) {
  # first, use svd to compute m-dimensional subspace
  # then, run naive()
  warning("Not yet implemented")
  nu_proj <- nu
  de_proj <- de
  return(naive(nu_proj, de_proj, ...))
}
