#' Least-squares heterodistributional subspace search
#'
#' @param df_numerator \code{data.frame} with exclusively numeric variables with
#' the numerator samples
#' @param df_denominator \code{data.frame} with exclusively numeric variables
#' with the denominator samples (must have the same variables as
#' \code{df_denominator})
#' @param m Scalar indicating the dimensionality of the reduced subspace
#' @param intercept \code{logical} Indicating whether to include an intercept
#' term in the model. Defaults to \code{TRUE}.
#' @param nsigma Integer indicating the number of sigma values (bandwidth
#' parameter of the Gaussian kernel gram matrix) to use in cross-validation.
#' @param sigma_quantile \code{NULL} or numeric vector with probabilities to
#' calculate the quantiles of the distance matrix to obtain sigma values. If
#' \code{NULL}, \code{nsigma} values between \code{0.25} and \code{0.75} are
#' used.
#' @param sigma \code{NULL} or a scalar value to determine the bandwidth of the
#' Gaussian kernel gram matrix. If \code{NULL}, \code{nsigma} values between
#' \code{0.25} and \code{0.75} are used.
#' @param nlambda Integer indicating the number of \code{lambda} values
#' (regularization parameter), by default, \code{lambda} is set to
#' \code{10^seq(3, -3, length.out = nlambda)}.
#' @param lambda \code{NULL} or numeric vector indicating the lambda values to
#' use in cross-validation
#' @param ncenters Maximum number of Gaussian centers in the kernel gram
#' matrix. Defaults to all numerator samples.
#' @param centers Numeric matrix with the same variables as \code{nu} and
#' \code{de} that are used as Gaussian centers in the kernel Gram matrix. By
#' default, the matrix \code{nu} is used as the matrix with Gaussian centers.
#' @param maxit Maximum number of iterations in the updating scheme.
#' @param parallel logical indicating whether to use parallel processing in the
#' cross-validation scheme.
#' @param nthreads \code{NULL} or integer indicating the number of threads to
#' use for parallel processing. If parallel processing is enabled, it defaults
#' to the number of available threads minus one.
#' @param progressbar Logical indicating whether or not to display a progressbar.
#' @export
#' @return \code{lhss}-object, containing all information to calculate the
#' density ratio using optimal sigma, optimal lambda and optimal weights.
#'
#' @return \code{lhss} returns \code{rhat}, the estimated density ratio.
#'
#' @examples
#' set.seed(1)
#' N <- 250
#' X <- cbind(rnorm(N), rnorm(N, 0, 0.5))
#' Y <- cbind(rnorm(N), sample(rep(c(-1, 1), times = N/2)) + rnorm(N))
#' out <- lhss(X, Y, m = 1, ncenters = 100)


lhss <- function(df_numerator, df_denominator, m = NULL, intercept = TRUE,
                     nsigma = 10, sigma_quantile = NULL, sigma = NULL,
                     nlambda = 10, lambda = NULL, ncenters = 200,
                     centers = NULL, maxit = 200, parallel = FALSE,
                     nthreads = NULL, progressbar = TRUE) {

  cl   <- match.call()
  nu   <- as.matrix(df_numerator)
  de   <- as.matrix(df_denominator)
  p    <- ncol(nu)
  n_nu <- nrow(nu)
  n_de <- nrow(de)

  check.dataform(nu, de)
  centers        <- check.centers(nu, centers, ncenters)
  symmetric      <- check.symmetric(nu, centers)
  parallel       <- check.parallel(parallel, nthreads, sigma, lambda)
  nthreads       <- check.threads(parallel, nthreads)
  sigma_quantile <- check.sigma_quantile.lhss(nsigma, sigma, sigma_quantile)
  sigma          <- if (is.null(sigma_quantile)) sigma else sigma_quantile
  is_quantile    <- !is.null(sigma_quantile)
  lambda         <- check.lambda(nlambda, lambda)
  m              <- check.subspace(m, p)

  res <- lhss_compute_alpha(nu, de, centers, symmetric, m, intercept, sigma,
                            is_quantile, lambda, maxit, parallel, nthreads,
                            progressbar)

  colnames(res$loocv) <- names(lambda) <- paste0("lambda", 1:length(lambda))
  rownames(res$loocv) <- names(sigma) <- paste0("sigma", 1:length(sigma))
  dimnames(res$sigmaopt) <- dimnames(res$loocv)
  dimnames(res$alpha) <- list(NULL, names(sigma), names(lambda))
  U <- array(res$Uopt, dim = c(p, m, length(sigma), length(lambda)),
             dimnames = list(NULL, NULL, names(sigma), names(lambda)))
  min_score <- arrayInd(which.min(res$loocv), dim(res$loocv))

  out <- list(
    df_numerator = df_numerator,
    df_denominator = df_denominator,
    alpha = res$alpha,
    cv_score = res$loocv,
    sigma = res$sigmaopt,
    sigma_quantiles = sigma_quantile,
    lambda = lambda,
    U = U,
    m = m,
    centers = centers,
    alpha_opt = res$alpha[, min_score[1], min_score[2]],
    lambda_opt = lambda[min_score[2]],
    sigma_opt = res$sigmaopt[min_score[1], min_score[2]],
    U_opt = U[,,min_score[1], min_score[2], drop = FALSE],
    call = cl
  )
  class(out) <- "lhss"
  out
}
