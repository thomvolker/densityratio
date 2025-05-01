#' Obtain predicted density ratio values from a \code{ulsif} object
#' @rdname predict.ulsif
#' @method predict ulsif
#' @param object A \code{ulsif} object
#' @param newdata Optional \code{matrix} new data set to compute the density
#' @param sigma A scalar with the Gaussian kernel width
#' @param lambda A scalar with the regularization parameter
#' @param ... Additional arguments to be passed to the function
#'
#' @return An array with predicted density ratio values from possibly new data,
#' but otherwise the numerator samples.
#'
#' @keywords predict ulsif
#' @seealso \code{\link{predict}}, \code{\link{ulsif}}
#'
#' @export
#'
#' @example inst/examples/ulsif-example.R

predict.ulsif <- function(object, newdata = NULL, sigma = c("sigmaopt", "all"), lambda = c("lambdaopt", "all"), ...) {
  newsigma <- check.sigma.predict(object, sigma)
  newlambda <- check.lambda.predict(object, lambda)
  newdata <- check.newdata(object, newdata)

  alpha <- extract_params(object, sigma = newsigma, lambda = newlambda, ...)
  nsigma <- length(newsigma)
  nlambda <- length(newlambda)
  dratio <- array(0, c(nrow(newdata), nsigma, nlambda))
  intercept <- nrow(object$alpha) > nrow(object$centers)

  for (i in 1:nsigma) {
    K <- kernel_gaussian(
      distance(newdata, object$model_matrices$ce, intercept),
      newsigma[i]
    )
    for (j in 1:nlambda) {
      dratio[, i, j] <- K %*% alpha[, i, j]
    }
  }
  dratio
}

#' Obtain predicted density ratio values from a \code{kliep} object
#' @rdname predict.kliep
#' @method predict kliep
#' @param object A \code{kliep} object
#' @param newdata Optional \code{matrix} new data set to compute the density
#' @param sigma A scalar with the Gaussian kernel width
#' @param ... Additional arguments to be passed to the function
#'
#' @return An array with predicted density ratio values from possibly new data,
#' but otherwise the numerator samples.
#'
#' @keywords predict kliep
#' @seealso \code{\link{predict}}, \code{\link{kliep}}
#'
#' @export
#' @example inst/examples/kliep-example.R


predict.kliep <- function(object, newdata = NULL, sigma = c("sigmaopt", "all"), ...) {
  newsigma <- check.sigma.predict(object, sigma)
  newdata <- check.newdata(object, newdata)

  alpha <- extract_params(object, sigma = newsigma, ...)
  nsigma <- ncol(alpha)
  dratio <- matrix(0, nrow(newdata), nsigma)
  intercept <- nrow(object$alpha) > nrow(object$centers)

  for (i in 1:nsigma) {
    K <- kernel_gaussian(
      distance(newdata, object$model_matrices$ce, intercept),
      newsigma[i]
    )
    dratio[, i] <- K %*% alpha[, i]
  }
  dratio
}

#' Obtain predicted density ratio values from a \code{kmm} object
#' @rdname predict.kmm
#' @method predict kmm
#' @param object A \code{kmm} object
#' @param newdata Optional \code{matrix} new data set to compute the density
#' @param sigma A scalar with the Gaussian kernel width
#' @param ... Additional arguments to be passed to the function
#'
#' @return An array with predicted density ratio values from possibly new data,
#' but otherwise the numerator samples.
#'
#' @keywords predict kmm
#' @seealso \code{\link{predict}}, \code{\link{kmm}}
#'
#' @export
#' @example inst/examples/kmm-example.R


predict.kmm <- function(object, newdata = NULL, sigma = c("sigmaopt", "all"), ...) {
  newsigma <- check.sigma.predict(object, sigma)
  newdata <- check.newdata(object, newdata)

  alpha <- extract_params(object, sigma = newsigma, ...)
  nsigma <- ncol(alpha)
  dratio <- matrix(0, nrow(newdata), nsigma)
  intercept <- nrow(object$alpha) > nrow(object$centers)

  for (i in 1:nsigma) {
    K <- kernel_gaussian(
      distance(newdata, object$model_matrices$ce, intercept),
      newsigma[i]
    )
    dratio[, i] <- K %*% alpha[, i]
  }
  dratio
}

#' Obtain predicted density ratio values from a \code{lhss} object
#' @rdname predict.lhss
#' @method predict lhss
#' @param object A \code{lhss} object
#' @param newdata Optional \code{matrix} new data set to compute the density
#' @param sigma A scalar with the Gaussian kernel width
#' @param lambda A scalar with the regularization parameter
#' @param ... Additional arguments to be passed to the function
#'
#' @return An array with predicted density ratio values from possibly new data,
#' but otherwise the numerator samples.
#'
#' @keywords predict lhss
#' @seealso \code{\link{predict}}, \code{\link{lhss}}
#'
#' @export
#' @example inst/examples/lhss-example.R


predict.lhss <- function(object, newdata = NULL, sigma = c("sigmaopt", "all"), lambda = c("lambdaopt", "all"), ...) {
  newlambda <- check.lambda.predict(object, lambda)
  lambdaind <- match(newlambda, object$lambda)
  lambdasigma <- check.lambdasigma.predict(object, sigma, newlambda, lambdaind)
  newdata <- check.newdata(object, newdata)

  alpha_U_sigma <- extract_params(object, lambda = newlambda, lambdasigma = lambdasigma, ...)
  nlambda <- length(newlambda)
  nsigma <- nrow(lambdasigma) / nlambda
  dratio <- array(0, c(nrow(newdata), nsigma, nlambda))
  intercept <- nrow(object$alpha) > nrow(object$centers)
  for (i in 1:nlambda) {
    for (j in 1:nsigma) {
      U_new <- alpha_U_sigma$U[, , j, i]
      alpha_new <- alpha_U_sigma$alpha[, j, i]
      K <- kernel_gaussian(
        distance(newdata %*% U_new, object$model_matrices$ce %*% U_new, intercept),
        sigma = alpha_U_sigma$sigma[j, i]
      )
      dratio[, j, i] <- K %*% alpha_new
    }
  }
  dratio
}

#' Obtain predicted density ratio values from a \code{spectral} object
#' @rdname predict.spectral
#' @method predict spectral
#' @param object A \code{spectral} object
#' @param newdata Optional \code{matrix} new data set to compute the density
#' @param sigma A scalar with the Gaussian kernel width
#' @param m integer indicating the dimension of the eigenvector expansion
#' @param ... Additional arguments to be passed to the function
#'
#' @return An array with predicted density ratio values from possibly new data,
#' but otherwise the numerator samples.
#'
#' @keywords predict spectral
#' @seealso \code{\link{predict}}, \code{\link{spectral}}
#'
#' @export

predict.spectral <- function(object, newdata = NULL, sigma = c("sigmaopt", "all"),
                             m = c("opt", "all"), ...) {
  newsigma <- check.sigma.predict(object, sigma)
  newM <- check.subspace.spectral.predict(object, m)
  newdata <- check.newdata(object, newdata)

  alpha_eigen <- extract_params(object, sigma = newsigma, m = newM, ...)
  dratio <- array(0, dim = c(nrow(newdata), length(newM), length(newsigma)))

  for (i in 1:length(newsigma)) {
    K <- kernel_gaussian(
      distance(newdata, object$model_matrices$ce, FALSE),
      newsigma[i]
    )
    D <- diag(length(alpha_eigen$Evals[, i]))
    diag(D) <- sqrt(nrow(object$model_matrices$ce)) / alpha_eigen$Evals[, i]
    phihatpred <- K %*% alpha_eigen$Evecs[, , i] %*% D
    for (m in 1:length(newM)) {
      dratio[, m, i] <- phihatpred[, seq_len(newM[m]), drop = FALSE] %*% alpha_eigen$alpha[seq_len(newM[m]), i, drop = FALSE]
    }
  }
  dratio
}

#' Obtain predicted density ratio values from a \code{naivedensityratio} object
#' @rdname predict.naivedensityratio
#' @method predict naivedensityratio
#' @param object A \code{naive} object
#' @param newdata Optional \code{matrix} new data set to compute the density
#' @param log A logical indicating whether to return the log of the density ratio
#' @param tol Minimal density value to avoid numerical issues
#' @param ... Additional arguments to be passed to the function
#'
#' @return An array with predicted density ratio values from possibly new data,
#' but otherwise the numerator samples.
#'
#' @keywords predict naive
#' @seealso \code{\link{predict}}, \code{\link{naive}}
#' @importFrom stats predict
#'
#' @export
#' @example inst/examples/naive-example.R


predict.naivedensityratio <- function(object, newdata = NULL, log = FALSE,
                                      tol = 1e-6, ...) {
  newdata <- check.newdata(object, newdata)
  newdata_proj <- asplit(predict(object$fit, newdata), 2)

  log_densities_nu <- mapply(
    function(density, newdata) {
      log(stats::approx(density, xout = newdata, yleft = tol, yright = tol)$y)
    }, object$density_numerator, newdata_proj
  )

  log_densities_de <- mapply(
    function(kde, newdata) {
      log(stats::approx(kde, xout = newdata, yleft = tol, yright = tol)$y)
    }, object$density_denominator, newdata_proj
  )

  densest_nu <- rowSums(log_densities_nu)
  densest_de <- rowSums(log_densities_de)

  if (log) {
    res <- densest_nu - densest_de
  } else {
    res <- exp(densest_nu - densest_de)
  }
  res
}


#' Extract parameters
#' @keywords internal

extract_params <- function(object, ...) {
  UseMethod("extract_params")
}

#' Obtain parameters from a \code{kliep} object
#'
#' @method extract_params kliep
#' @keywords internal

extract_params.kliep <- function(object, sigma, ...) {
  if (all(sigma %in% object$sigma)) {
    which_sigma <- match(sigma, object$sigma)
    alpha <- object$alpha[, which_sigma, drop = FALSE]
  } else {
    alpha <- update(object, sigma = sigma, cv = FALSE, ...)$alpha
  }
  alpha
}

#' Obtain parameters from a \code{kmm} object
#'
#' @method extract_params kmm
#' @keywords internal

extract_params.kmm <- extract_params.kliep

#' Obtain parameters from a \code{ulsif} object
#'
#' @method extract_params ulsif
#' @keywords internal

extract_params.ulsif <- function(object, sigma, lambda, ...) {
  if (all(sigma %in% object$sigma) & all(lambda %in% object$lambda)) {
    which_sigma <- match(sigma, object$sigma)
    which_lambda <- match(lambda, object$lambda)
    alpha <- object$alpha[, which_sigma, which_lambda, drop = FALSE]
  } else {
    alpha <- update(object, sigma = sigma, lambda = lambda, ...)$alpha
  }
  alpha
}

#' Obtain parameters from a \code{lhss} object
#'
#' @method extract_params lhss
#' @keywords internal

extract_params.lhss <- function(object, lambda, lambdasigma, ...) {
  alpha_old <- object$alpha
  nlambda <- length(lambda)
  nsigma <- nrow(lambdasigma) / nlambda
  alpha <- array(0, dim = c(dim(alpha_old)[1], nsigma, nlambda))
  U <- array(0, dim = c(nrow(object$U_opt), ncol(object$U_opt), nsigma, nlambda))
  sigma <- matrix(0, nsigma, nlambda)
  for (i in 1:nlambda) {
    for (j in 1:nsigma) {
      ind <- (i - 1) * nsigma + j
      if (!is.na(lambdasigma[ind, 1]) & !is.na(lambdasigma[ind, 2])) {
        alpha[, j, i] <- object$alpha[, lambdasigma[ind, 2], lambdasigma[ind, 1]]
        U[, , j, i] <- object$U[, , lambdasigma[ind, 2], lambdasigma[ind, 1]]
        sigma[j, i] <- object$sigma[lambdasigma[ind, 2], lambdasigma[ind, 1]]
      } else {
        if (is.na(lambdasigma[ind, 4])) {
          fitnew <- update(
            object,
            sigma = NULL,
            sigma_quantile = object$sigma_quantiles[j],
            lambda = lambdasigma[ind, 3],
            ...
          )
        } else {
          fitnew <- update(
            object,
            sigma = lambdasigma[ind, 4],
            lambda = lambdasigma[ind, 3],
            ...
          )
        }
        alpha[, j, i] <- fitnew$alpha
        U[, , j, i] <- fitnew$U_opt
        sigma[j, i] <- fitnew$sigma_opt
      }
    }
  }
  list(alpha = alpha, U = U, sigma = sigma)
}

#' Obtain parameters from a \code{spectral} object
#'
#' @method extract_params spectral
#' @keywords internal

extract_params.spectral <- function(object, sigma, m, ...) {
  maxM <- max(m) # largest subspace dimension
  nsigma <- length(sigma) # number of new sigma values in predict
  which_sigma <- match(sigma, object$sigma) # indices of original sigma values per new sigma

  if (all(sigma %in% object$sigma) & maxM <= max(object$m)) {
    # extract parameters
    Evals <- object$Evals[1:maxM, which_sigma, drop = FALSE]
    Evecs <- object$Evecs[, 1:maxM, which_sigma, drop = FALSE]
    alpha <- object$alpha[1:maxM, which_sigma, drop = FALSE]
  } else {
    if (maxM <= max(object$m)) {
      sigma_new <- sigma[which(is.na(which_sigma))]
      sigma_old_ind <- which_sigma[!is.na(which_sigma)]

      # update fit object to accomodate new parameters
      fit <- update(object, sigma = sigma_new, m = maxM, cv = FALSE, ...)

      # initialize empty objects to store results
      alpha <- matrix(0, maxM, nsigma)
      Evals <- matrix(0, maxM, nsigma)
      Evecs <- array(0, dim = c(dim(object$Evecs)[1], maxM, nsigma))

      # store old results on correct places
      alpha[, sigma_old_ind] <- object$alpha[1:maxM, sigma_old_ind]
      Evals[, sigma_old_ind] <- object$Evals[1:maxM, sigma_old_ind]
      Evecs[, , sigma_old_ind] <- object$Evecs[, 1:maxM, sigma_old_ind]

      # store new results on correct places
      alpha[, which(is.na(which_sigma))] <- fit$alpha
      Evals[, which(is.na(which_sigma))] <- fit$Evals
      Evecs[, , which(is.na(which_sigma))] <- fit$Evecs
    } else {
      fit <- update(object, sigma = sigma, m = m, cv = FALSE, ...)
      alpha <- fit$alpha
      Evals <- fit$Evals
      Evecs <- fit$Evecs
    }
  }
  list(alpha = alpha, Evals = Evals, Evecs = Evecs)
}
