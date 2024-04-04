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
#' @examples
#' x <- rnorm(100) |> matrix(100)
#' y <- rnorm(200, 1, 2) |> matrix(200)
#' fit1 <- ulsif(x, y)
#' predict(fit1)
#' predict(fit1, newdata = rbind(x, y))
#' predict(fit1, newdata = rbind(x, y), sigma = 2, lambda = 3)

predict.ulsif <- function(object, newdata = NULL, sigma = c("sigmaopt", "all"), lambda = c("lambdaopt", "all"), ...) {

  newsigma  <- check.sigma.predict(object, sigma)
  newlambda <- check.lambda.predict(object, lambda)
  newdata <- check.newdata(object, newdata)

  alpha     <- extract_params(object, sigma = newsigma, lambda = newlambda, ...)
  nsigma    <- length(newsigma)
  nlambda   <- length(newlambda)
  dratio    <- array(0, c(nrow(newdata), nsigma, nlambda))
  intercept <- nrow(object$alpha) > nrow(object$centers)

  for (i in 1:nsigma) {
    K <- distance(newdata, object$centers, intercept) |> kernel_gaussian(newsigma[i])
    for (j in 1:nlambda) {
      dratio[ , i, j] <- K %*% alpha[, i, j]
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
#'
#' @examples
#' x <- rnorm(100) |> matrix(100)
#' y <- rnorm(200, 1, 2) |> matrix(200)
#' fit1 <- kliep(x, y)
#' predict(fit1)
#' predict(fit1, newdata = rbind(x, y))
#' predict(fit1, newdata = rbind(x, y), sigma = 2)

predict.kliep <- function(object, newdata = NULL, sigma = c("sigmaopt", "all"), ...) {

  newsigma <- check.sigma.predict(object, sigma)
  newdata  <- check.newdata(object, newdata)

  alpha <- extract_params(object, sigma = newsigma, ...)
  nsigma <- ncol(alpha)
  dratio <- matrix(0, nrow(newdata), nsigma)
  intercept <- nrow(object$alpha) > nrow(object$centers)

  for (i in 1:nsigma) {
    K <- distance(newdata, object$centers, intercept) |> kernel_gaussian(newsigma[i])
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
#'
#' @examples
#' x <- rnorm(100) |> matrix(50)
#' y <- rnorm(200, 1, 2) |> matrix(100)
#' fit1 <- lhss(x, y, m = 1)
#' predict(fit1)
#' predict(fit1, newdata = rbind(x, y))
#' predict(fit1, newdata = rbind(x, y), sigma = 2)

predict.lhss <- function(object, newdata = NULL, sigma = c("sigmaopt", "all"), lambda = c("lambdaopt", "all"), ...) {

  newlambda   <- check.lambda.predict(object, lambda)
  lambdaind   <- match(newlambda, object$lambda)
  lambdasigma <- check.lambdasigma.predict(object, sigma, newlambda, lambdaind)
  newdata     <- check.newdata(object, newdata)

  alpha_U_sigma <- extract_params(object, lambda = newlambda, lambdasigma = lambdasigma, ...)
  nlambda     <- length(newlambda)
  nsigma      <- nrow(lambdasigma) / nlambda
  dratio    <- array(0, c(nrow(newdata), nsigma, nlambda))
  intercept <- nrow(object$alpha) > nrow(object$centers)
  for (i in 1:nlambda) {
    for (j in 1:nsigma) {
      U_new <- alpha_U_sigma$U[ , , j, i]
      alpha_new <- alpha_U_sigma$alpha[ , j, i]
      K <- distance(newdata %*% U_new, object$centers %*% U_new, intercept) |>
        kernel_gaussian(sigma = alpha_U_sigma$sigma[j, i])
      dratio[ , j, i] <- K %*% alpha_new
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
#' @param J integer indicating the dimension of the eigenvector expansion
#' @param tol A scalar indicating the smallest eligible density ratio value
#' (used to censor negative predicted density ratio values).
#' @param ... Additional arguments to be passed to the function
#'
#' @return An array with predicted density ratio values from possibly new data,
#' but otherwise the numerator samples.
#'
#' @keywords predict spectral
#' @seealso \code{\link{predict}}, \code{\link{spectral}}
#'
#' @export
#'
#' @examples
#' x <- rnorm(100) |> matrix(100)
#' y <- rnorm(200, 1, 2) |> matrix(200)
#' fit1 <- spectral(x, y)
#' predict(fit1)
#' predict(fit1, newdata = rbind(x, y))
#' predict(fit1, newdata = rbind(x, y), sigma = 2, J = 10)

predict.spectral <- function(object, newdata = NULL, sigma = c("sigmaopt", "all"),
                             J = c("Jopt", "all"), tol = 1e-6, ...) {

  newsigma <- check.sigma.predict(object, sigma)
  newJ     <- check.J.predict(object, J)
  newdata  <- check.newdata(object, newdata)

  alpha_eigen <- extract_params(object, sigma = newsigma, J = newJ, ...)
  dratio <- array(0, dim = c(nrow(newdata), length(newJ), length(newsigma)))

  for (i in 1:length(newsigma)) {
    K <- distance(newdata, object$centers, FALSE) |> kernel_gaussian(newsigma[i])
    phihatpred <- K %*% alpha_eigen$Evecs[,,i] %*% diag(sqrt(nrow(object$centers))/alpha_eigen$Evals[,i])
    for (j in 1:length(newJ)) {
      dratio[ , j, i] <- pmax(tol, phihatpred[,seq_len(newJ[j])] %*% alpha_eigen$alpha[seq_len(newJ[j]),i])
    }
  }
  dratio
}

#' Predict function for density object
#'
#' @method predict density
#' @keywords internal
#' @importFrom stats approx

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
#' @rdname predict.naivedensityratio
#' @method predict naivedensityratio
#' @param object A \code{naive} object
#' @param newdata Optional \code{matrix} new data set to compute the density
#' @param log A logical indicating whether to return the log of the density ratio
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
#'
#' @examples
#' x <- rnorm(100) |> matrix(100)
#' y <- rnorm(200, 1, 2) |> matrix(200)
#' fit1 <- naive(x, y)
#' predict(fit1)
#' predict(fit1, newdata = rbind(x, y))
#' predict(fit1, newdata = rbind(x, y), log = TRUE)

predict.naivedensityratio <- function(object, newdata = NULL, log = FALSE, ...) {
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

  # return log-density difference (or its exponent, i.e., the density ratio)
  ld_dif <- ld_nu - ld_de
  if (log)  {
    res <- ld_dif
  } else {
    res <- exp(ld_dif)
  }
  res
}

#' Obtain predicted density ratio values from a \code{naivesubspace} object
#' @rdname predict.naivesubspacedensityratio
#' @method predict naivesubspacedensityratio
#' @param object A \code{naivesubspace} object
#' @param newdata Optional \code{matrix} new data set to compute the density
#' @param log A logical indicating whether to return the log of the density ratio
#' @param ... Additional arguments to be passed to the function
#'
#' @return An array with predicted density ratio values from possibly new data,
#' but otherwise the numerator samples.
#'
#' @keywords predict naivedensityratio
#' @seealso \code{\link{predict}}, \code{\link{naivesubspace}}
#'
#' @export
#'
#' @examples
#' x <- rnorm(100) |> matrix(100)
#' y <- rnorm(200, 1, 2) |> matrix(200)
#' fit1 <- naivesubspace(x, y)
#' predict(fit1)
#' predict(fit1, newdata = rbind(x, y))
#' predict(fit1, newdata = rbind(x, y), log = TRUE)

predict.naivesubspacedensityratio <- function(object, newdata = NULL, log = FALSE, ...) {

  newdata <- check.newdata(object, newdata)

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

  # return log-density difference (or its exponent, i.e., the density ratio)
  ld_dif <- ld_nu - ld_de
  if (log)  {
    res <- ld_dif
  } else {
    res <- exp(ld_dif)
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
    alpha <- object$alpha[ , which_sigma, drop = FALSE]
  } else {
    alpha <- update(object, sigma = sigma, cv = FALSE, ...)$alpha
  }
  alpha
}

#' Obtain parameters from a \code{ulsif} object
#'
#' @method extract_params ulsif
#' @keywords internal

extract_params.ulsif <- function(object, sigma, lambda, ...) {
  if (all(sigma %in% object$sigma) & all(lambda %in% object$lambda)) {
    which_sigma <- match(sigma, object$sigma)
    which_lambda <- match(lambda, object$lambda)
    alpha <- object$alpha[ , which_sigma, which_lambda, drop = FALSE]
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
      ind <- (i-1) * nsigma + j
      if (!is.na(lambdasigma[ind, 1]) & !is.na(lambdasigma[ind, 2])) {
        alpha[ , j, i] <- object$alpha[, lambdasigma[ind, 2], lambdasigma[ind, 1]]
        U[ , , j, i] <- object$U[, , lambdasigma[ind, 2], lambdasigma[ind, 1]]
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
        alpha[ , j, i] <- fitnew$alpha
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

extract_params.spectral <- function(object, sigma, J, ...) {

  maxJ <- max(J) #largest subspace dimension
  nsigma <- length(sigma) # number of new sigma values in predict
  which_sigma <- match(sigma, object$sigma) # indices of original sigma values per new sigma

  if (all(sigma %in% object$sigma) & all(J <= max(object$J))) {
    # extract parameters
    Evals <- object$Evals[1:maxJ, which_sigma, drop = FALSE]
    Evecs <- object$Evecs[ , 1:maxJ, which_sigma, drop = FALSE]
    alpha <- object$alpha[1:maxJ, which_sigma, drop = FALSE]
  } else {
    if (maxJ <= max(object$J)) {
      sigma_new <- sigma[which(is.na(which_sigma))]
      sigma_old_ind <- which_sigma[!is.na(which_sigma)]

      # update fit object to accomodate new parameters
      fit <- update(object, sigma = sigma_new, J = maxJ, cv = FALSE, ...)

      # initialize empty objects to store results
      alpha <- matrix(0, maxJ, nsigma)
      Evals <- matrix(0, maxJ, nsigma)
      Evecs <- array(0, dim = c(dim(object$Evecs)[1], maxJ, nsigma))

      # store old results on correct places
      alpha[, sigma_old_ind]   <- object$alpha[1:max(J), sigma_old_ind]
      Evals[, sigma_old_ind]   <- object$Evals[1:maxJ, sigma_old_ind]
      Evecs[, , sigma_old_ind] <- object$Evecs[, 1:max(J), sigma_old_ind]

      # store new results on correct places
      alpha[, which(is.na(which_sigma))] <- fit$alpha
      Evals[, which(is.na(which_sigma))] <- fit$Evals
      Evecs[, , which(is.na(which_sigma))] <- fit$Evecs
    } else {
      fit <- update(object, sigma = sigma, J = J, cv = FALSE, ...)
      alpha <- fit$alpha
      Evals <- fit$Evals
      Evecs <- fit$Evecs
    }
  }
  list(alpha = alpha, Evals = Evals, Evecs = Evecs)
}

