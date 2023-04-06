#' Predict function for density ratio estimation
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

#' Naive density ratio estimation
#'
#' @param nu Numeric matrix with numerator samples
#' @param de Numeric matrix with denominator samples (must have the same
#' variables as \code{nu})
#' @param ... arguments passed to [stats::density]
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
  # project de to m-dimensional space
  de_centered <- scale(de, scale = FALSE)
  V <- svd(de_centered, nu = m, nv = m)$v
  de_proj <- de_centered %*% V

  nu_centered <- scale(nu, center = attr(de_centered, "scaled:center"), scale = FALSE)
  nu_proj <- nu_centered %*% V

  return(naive(nu_proj, de_proj, ...))
}
