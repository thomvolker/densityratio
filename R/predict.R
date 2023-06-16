
predict.ulsif <- function(object, newdata = NULL, sigma = c("sigmaopt", "all"), lambda = c("lambdaopt", "all")) {

  newsigma  <- check.sigma.predict(object, sigma)
  newlambda <- check.lambda.predict(object, lambda)
  newdata <- check.newdata(object, newdata)

  alpha   <- extract.alpha(object, newsigma, newlambda)
  nlambda <- dim(alpha)[2]
  nsigma  <- dim(alpha)[3]
  dratio  <- array(0, c(nrow(alpha), nlambda, nsigma))

  for (i in 1:nsigma) {
    K <- distance(newdata, object$centers) |> kernel_gaussian(newsigma[i])
    for (j in 1:length(lambda)) {
      dratio[ , j, i] <- K %*% alpha[, j, i]
    }
  }
  dratio
}

predict.kliep <- function(object, newdata = NULL, sigma = c("sigmaopt", "all")) {
  newsigma <- check.sigma.predict(object, sigma)
  newdata  <- check.newdata(object, newdata)

  alpha <- extract.alpha(object, newsigma, lambda = NULL)
  nsigma <- ncol(alpha)
  dratio <- matrix(0, nrow(alpha), nsigma)

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
      alpha <- ulsif(object$numerator, object$denominator, sigma = sigma,
                     lambda = lambda, centers = object$centers, progressbar = FALSE)$alpha
    }
  }
  alpha
}



