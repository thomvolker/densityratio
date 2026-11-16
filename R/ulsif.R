#' Relative unconstrained least-squares importance fitting (rulsif)
#' 
#' Estimate the a-relative density ratio 
#' \deqn{r_a(x) = \frac{p_{nu}(x)}{a p_{nu}(x) + (1-a) p_{de}(x)}} with 
#' relative ulsif. \code{ulsif()} is a special case of \code{rulsif()} with
#' \code{a = 0}, which estimates the ordinary density ratio
#' \eqn{p_{nu}(x) / p_{de}(x)} and is kept for backward compatibility: 
#' \code{ulsif(...)} is equivalent to \code{rulsif(..., a = 0)}. 
#' 
#' @details
#' For \code{a > 0} the relative density ratio is bounded above by \code{1/a},
#' which makes the estimator more stable when the denominator density is small
#' where the numerator density is not (Yamada et al., 2013). 
#' 
#' The bandwidth \code{sigma} and regularization parameter \code{lambda} are 
#' chosen by leave-one-out cross-validation, computed in closed-form. Leaving
#' out one numerator and one denominator sample changes the system by a rank-2
#' 2, which is handled with the Sherman-Morrison-Woodbury identity. For 
#' \code{a = 0}, this reduces to the rank-1 update of Kanamori et al. (2009). 
#'
#' @param df_numerator \code{data.frame} with exclusively numeric variables with
#' the numerator samples
#' @param df_denominator \code{data.frame} with exclusively numeric variables
#' with the denominator samples (must have the same variables as
#' \code{df_numerator})
#' @param a Numeric scaler with \code{0 <= a < 1}: the mixing proportion of
#' the numerator in the reference density \code{a p_{nu} + (1-a)p_{de}}. 
#' \code{a = 0} gives ordinary ulsif. Defaults to \code{a = 0.1} in 
#' \code{rulsif()}.
#' @param intercept \code{logical} Indicating whether to include an intercept
#' term in the model. Defaults to \code{FALSE}.
#' @param scale \code{"numerator"}, \code{"denominator"}, or \code{NULL},
#' indicating whether to standardize each numeric variable according to the
#' numerator means and standard deviations, the denominator means and standard
#' deviations, or apply no standardization at all.
#' @param nsigma Integer indicating the number of sigma values (bandwidth
#' parameter of the Gaussian kernel gram matrix) to use in cross-validation.
#' @param sigma_quantile \code{NULL} or numeric vector with probabilities to
#' calculate the quantiles of the distance matrix to obtain sigma values. If
#' \code{NULL}, \code{nsigma} values between \code{0.05} and \code{0.95} are
#' used.
#' @param sigma \code{NULL} or a scalar value to determine the bandwidth of the
#' Gaussian kernel gram matrix. If \code{NULL}, \code{nsigma} values between
#' \code{0.05} and \code{0.95} are used.
#' @param nlambda Integer indicating the number of \code{lambda} values
#' (regularization parameter), by default, \code{lambda} is set to
#' \code{10^seq(3, -3, length.out = nlambda)}.
#' @param lambda \code{NULL} or numeric vector indicating the lambda values to
#' use in cross-validation
#' @param ncenters \code{integer} Maximum number of Gaussian centers in the kernel 
#' gram matrix.
#' @param centers \code{NULL} or \code{data.frame} with the same dimensions as the
#' data, indicating the centers for the Gaussian kernel gram matrix.
#' @param parallel logical indicating whether to use parallel processing in the
#' cross-validation scheme.
#' @param nthreads \code{NULL} or integer indicating the number of threads to
#' use for parallel processing. If parallel processing is enabled, it defaults
#' to the number of available threads minus one.
#' @param progressbar Logical indicating whether or not to display a progressbar.
#' @references Kanamori, T., Hido, S., & Sugiyama, M. (2009). A least-squares
#' approach to direct importance estimation. Journal of Machine Learning
#' Research, 10, 1391-1445. Available from
#' <https://jmlr.org/papers/v10/kanamori09a.html>
#' Yamada, M., Suzuki, T., Kanamori, T., Hachiya, H., & Sugiyama, M. (2013). 
#' Relative density-ratio estimation for robust distribution comparison. 
#' Neural Computation, 25, 1324-1370. <doi:10.1162/NECO_a_00442>
#' @seealso \code{\link{summary.ulsif}}, \code{\link{predict.ulsif}}.
#'
#' @return An object of class \code{c("rulsif", "ulsif")} from \code{rulsif()},
#' or of class \code{"ulsif"} from \code{ulsif()}, containing the fitted model, 
#' the cross-validation scores and the optimal \code{sigma} and \code{lambda}.
#' The element \code{a} stores the mixing proportion. 
#'
#' @example inst/examples/ulsif-example.R
#' @rdname ulsif
#' @export



rulsif <- function(df_numerator, df_denominator, a = 0.1, intercept = FALSE, scale = "numerator",
                   nsigma = 10, sigma_quantile = NULL, sigma = NULL, nlambda = 20,
                   lambda = NULL, ncenters = 200, centers = NULL, parallel = FALSE, 
                   nthreads = NULL, progressbar = TRUE) {
 
  cl <- match.call()

  nu <- check.datatype(df_numerator)
  de <- check.datatype(df_denominator)

  check.variables(nu, de, centers)

  df_centers <- check.centers(rbind(nu, de), centers, ncenters)
  dat <- check.dataform(nu, de, df_centers, is.null(centers), NULL, scale)

  intercept <- check.intercept(intercept)

  check.a.rulsif(a)

  dist_nu <- distance(dat$nu, dat$ce, intercept)
  dist_de <- distance(dat$de, dat$ce, intercept)

  sigma <- check.sigma(nsigma, sigma_quantile, sigma, dist_nu)
  lambda <- check.lambda(nlambda, lambda)

  parallel <- check.parallel(parallel, nthreads, lambda)
  nthreads <- check.threads(parallel, nthreads)

  res <- compute_ulsif(dist_nu, dist_de, intercept, sigma, lambda, parallel, nthreads, progressbar, a = a)
  loocv_scores <- res$loocv_score
  colnames(loocv_scores) <- names(lambda) <- paste0("lambda", 1:length(lambda))
  rownames(loocv_scores) <- names(sigma) <- paste0("sigma", 1:length(sigma))
  dimnames(res$alpha) <- list(NULL, names(sigma), names(lambda))
  min_score <- arrayInd(which.min(loocv_scores), dim(loocv_scores))

  out <- list(
    df_numerator = df_numerator,
    df_denominator = df_denominator,
    a = a,
    alpha = res$alpha,
    cv_score = loocv_scores,
    scale = scale,
    sigma = sigma,
    lambda = lambda,
    centers = df_centers,
    model_matrices = dat,
    alpha_opt = res$alpha[, min_score[1], min_score[2]],
    lambda_opt = lambda[min_score[2]],
    sigma_opt = sigma[min_score[1]],
    call = cl
  )

  class(out) <- c("rulsif", "ulsif")
  out
}

#' @rdname ulsif
#' @export

ulsif <- function(df_numerator, df_denominator, intercept = FALSE, scale = "numerator",
                   nsigma = 10, sigma_quantile = NULL, sigma = NULL, nlambda = 20,
                   lambda = NULL, ncenters = 200, centers = NULL, parallel = FALSE, 
                   nthreads = NULL, progressbar = TRUE) {
  cl <- match.call()
  out <- rulsif(
    df_numerator = df_numerator, 
    df_denominator = df_denominator, 
    a = 0, 
    intercept = intercept, 
    scale = scale, 
    nsigma = nsigma, 
    sigma_quantile = sigma_quantile,
    sigma = sigma, 
    nlambda = nlambda, 
    lambda = lambda, 
    ncenters = ncenters, 
    centers = centers, 
    parallel = parallel, 
    nthreads = nthreads, 
    progressbar = progressbar
  )
  out$call <- cl
  class(out) <- "ulsif"
  out
}