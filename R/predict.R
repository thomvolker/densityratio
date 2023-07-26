#' Obtain predicted density ratio values from a \code{ulsif} object
#'
#' @rdname predict
#' @param object Object of class \code{ulsif} or \code{kliep}
#' @param newdata Optional \code{matrix} new data set to compute the density
#' ratio values of.
#' @param sigma Either one of c("sigmaopt", "all") for the optimal sigma value
#' or all sigma values used to far, or a scalar or numeric vector with new
#' sigma values to use.
#' @param lambda Either one of c("lambdaopt", "all") for the optimal lambda
#' value or all lambda values used so far, or a scalar or numeric vector with
#' new lambda values to use.
#' @param ... further arguments passed to or from other methods.
#' @return An array with predicted density ratio values from possibly new data,
#' but otherwise the numerator samples.
#' @method predict ulsif
#' @export


predict.ulsif <- function(object, newdata = NULL, sigma = c("sigmaopt", "all"), lambda = c("lambdaopt", "all"), ...) {

  newsigma  <- check.sigma.predict(object, sigma)
  newlambda <- check.lambda.predict(object, lambda)
  newdata <- check.newdata(object, newdata)

  alpha   <- extract.alpha(object, newsigma, newlambda)
  nlambda <- dim(alpha)[2]
  nsigma  <- dim(alpha)[3]
  dratio  <- array(0, c(nrow(newdata), nlambda, nsigma))

  for (i in 1:nsigma) {
    K <- distance(newdata, object$centers) |> kernel_gaussian(newsigma[i])
    for (j in 1:nlambda) {
      dratio[ , j, i] <- K %*% alpha[, j, i]
    }
  }
  dratio
}

#' Obtain predicted density ratio values from a \code{kliep} object
#'
#' @rdname predict
#' @method predict kliep
#' @export

predict.kliep <- function(object, newdata = NULL, sigma = c("sigmaopt", "all"), ...) {
  newsigma <- check.sigma.predict(object, sigma)
  newdata  <- check.newdata(object, newdata)

  alpha <- extract.alpha(object, newsigma, lambda = NULL)
  nsigma <- ncol(alpha)
  dratio <- matrix(0, nrow(newdata), nsigma)

  for (i in 1:nsigma) {
    K <- distance(newdata, object$centers) |> kernel_gaussian(newsigma[i])
    dratio[, i] <- K %*% alpha[, i]
  }
  dratio
}

extract.alpha <- function(object, sigma, lambda) {

  if (inherits(object, "kliep")) {
    if (all(sigma %in% object$sigma)) {
      which_sigma <- which(object$sigma %in% sigma)
      alpha <- object$alpha[ , which_sigma, drop = FALSE]
    } else {
      alpha <- kliep(object$df_numerator, object$df_denominator, sigma = sigma,
                     centers = object$centers, cv = FALSE, epsilon = object$epsilon,
                     maxit = object$maxit, progressbar = FALSE)$alpha
    }
  } else if (inherits(object, "ulsif")) {
    if (all(sigma %in% object$sigma) & all(lambda %in% object$lambda)) {
      which_sigma <- which(object$sigma %in% sigma)
      which_lambda <- which(object$lambda %in% lambda)
      alpha <- object$alpha[ , which_lambda, which_sigma, drop = FALSE]
    } else {
      alpha <- ulsif(object$df_numerator, object$df_denominator, sigma = sigma,
                     lambda = lambda, centers = object$centers, progressbar = FALSE)$alpha
    }
  }
  alpha
}


#' Predict function for density object
#'
#' @method predict density
#' @importFrom stats predict
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

#' Obtain predicted density ratio values from a \code{naivedensityratio} object
#'
#' @rdname predict
#' @method predict naivedensityratio
#' @importFrom stats approx predict
#'
#' @export
predict.naivedensityratio <- function(object, newdata = NULL, ...) {
  newdata <- check.newdata(object, newdata)
  P <- ncol(newdata)
  N <- nrow(newdata)

  # work on log scale
  # log-densities
  ld_nu <- numeric(N)
  ld_de <- numeric(N)

  # just add log-densities together
  for (p in 1:P) {
    ld_nu <- ld_nu + predict(object$density_numerator[[p]], newdata[,p], log = TRUE)
    ld_de <- ld_de + predict(object$density_denominator[[p]], newdata[,p], log = TRUE)
  }

  # return the exponent of the log-density difference, i.e., the density ratio
  return(exp(ld_nu - ld_de))
}

#' Obtain predicted density ratio values from a \code{naivesubspacedensityratio} object
#'
#' @rdname predict
#' @method predict naivesubspacedensityratio
#' @importFrom stats approx predict
#'
#' @export
predict.naivesubspacedensityratio <- function(object, newdata = NULL, ...) {
  newdata <- check.newdata(object, newdata)
  P <- ncol(newdata)
  N <- nrow(newdata)

  # project newdata to subspace
  nd <- as.matrix(newdata)
  nd_centered <- scale(nd, center = object$center, scale = FALSE)
  nd_proj <- nd_centered %*% object$projection_matrix

  # work on log scale
  # log-densities
  ld_nu <- numeric(N)
  ld_de <- numeric(N)

  # just add log-densities together
  for (p in 1:object$subspace_dim) {
    ld_nu <- ld_nu + predict(object$density_numerator[[p]], nd_proj[,p], log = TRUE)
    ld_de <- ld_de + predict(object$density_denominator[[p]], nd_proj[,p], log = TRUE)
  }

  # return the exponent of the log-density difference, i.e., the density ratio
  return(exp(ld_nu - ld_de))
}



