#' Kernel mean matching approach to density ratio estimation
#'
#' @param df_numerator \code{data.frame} with exclusively numeric variables with
#' the numerator samples
#' @param df_denominator \code{data.frame} with exclusively numeric variables
#' with the denominator samples (must have the same variables as
#' \code{df_denominator})
#' @param scale \code{"numerator"}, \code{"denominator"}, or \code{NULL},
#' indicating whether to standardize each numeric variable according to the
#' numerator means and standard deviations, the denominator means and standard
#' deviations, or apply no standardization at all.
#' @param constrained \code{logical} equals \code{FALSE} to use unconstrained
#' optimization, \code{TRUE} to use constrained optimization. Defaults to
#' \code{FALSE}.
#' @param nsigma Integer indicating the number of sigma values (bandwidth
#' parameter of the Gaussian kernel gram matrix) to use in cross-validation.
#' @param sigma_quantile \code{NULL} or numeric vector with probabilities to
#' calculate the quantiles of the distance matrix to obtain sigma values. If
#' \code{NULL}, \code{nsigma} values between \code{0.25} and \code{0.75} are
#' used.
#' @param sigma \code{NULL} or a scalar value to determine the bandwidth of the
#' Gaussian kernel gram matrix. If \code{NULL}, \code{nsigma} values between
#' \code{0.25} and \code{0.75} are used.
#' @param ncenters Maximum number of Gaussian centers in the kernel gram
#' matrix. Defaults to all numerator samples.
#' @param centers Option to specify the Gaussian samples manually.
#' @param cv Logical indicating whether or not to do cross-validation
#' @param nfold Number of cross-validation folds used in order to calculate the
#' optimal \code{sigma} value (default is 5-fold cv).
#' @param parallel logical indicating whether to use parallel processing in the
#' cross-validation scheme.
#' @param nthreads \code{NULL} or integer indicating the number of threads to
#' use for parallel processing. If parallel processing is enabled, it defaults
#' to the number of available threads minus one.
#' @param progressbar Logical indicating whether or not to display a progressbar.
#' @param osqp_settings Optional: settings to pass to the \code{osqp} solver for
#' constrained optimization.
#' @param cluster Optional: a cluster object to use for parallel processing,
#' see \code{parallel::makeCluster}.
#' @importFrom stats predict
#' @importFrom pbapply pblapply
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom osqp solve_osqp osqpSettings

#' @export
#'
#' @return \code{kmm}-object, containing all information to calculate the
#' density ratio using optimal sigma and optimal weights.
#'
#' @example inst/examples/kmm-example.R

kmm <- function(df_numerator, df_denominator, scale = "numerator",
                constrained = FALSE, nsigma = 10, sigma_quantile = NULL,
                sigma = NULL, ncenters = 200, centers = NULL, cv = TRUE,
                nfold = 5, parallel = FALSE, nthreads = NULL,
                progressbar = TRUE, osqp_settings = NULL, cluster = NULL) {
  cl <- match.call()
  nu <- check.datatype(df_numerator)
  de <- check.datatype(df_denominator)

  check.variables(nu, de)

  df_centers <- check.centers(rbind(nu, de), centers, ncenters)
  dat <- check.dataform(nu, de, df_centers, is.null(centers), NULL, scale)

  Dd <- distance(dat$de, dat$ce, FALSE)
  sigma <- check.sigma(nsigma, sigma_quantile, sigma, Dd)

  parallel <- check.parallel(parallel, nthreads, sigma)
  if (parallel && constrained) {
    if (is.null(cluster) && is.null(nthreads)) {
      nthreads <- parallel::detectCores() - 1
      cluster <- parallel::makeCluster(nthreads)
    } else if (is.null(cluster) && !is.null(nthreads)) {
      cluster <- parallel::makeCluster(nthreads)
    }
    on.exit(parallel::stopCluster(cluster))
  }

  cv_ind_nu <- check.nfold(cv, nfold, nrow(dat$nu))
  cv_ind_de <- check.nfold(cv, nfold, nrow(dat$de))

  if (is.null(osqp_settings)) {
    osqp_settings <- osqp::osqpSettings(verbose = FALSE)
  }
  if (constrained) {
    if (!progressbar) {
      pbapply::pboptions(type = "none")
    }
    constrained_out <- pbapply::pblapply(sigma, function(s) {
      compute_kmm(dat$nu, dat$de, dat$ce, Dd, s, cv_ind_nu, cv_ind_de,
        parallel = FALSE, nthreads = 0, progressbar = FALSE,
        constrained = TRUE, settings = osqp_settings
      )
    }, cl = cluster)
    res <- list(
      alpha = do.call(cbind, lapply(constrained_out, function(x) x$alpha)),
      loss = sapply(constrained_out, function(x) x$loss)
    )
  } else {
    nthreads <- check.threads(parallel, nthreads)
    res <- compute_kmm(
      dat$nu, dat$de, dat$ce, Dd, sigma, cv_ind_nu, cv_ind_de,
      parallel, nthreads, progressbar, constrained, osqp_settings
    )
  }

  out <- list(
    df_numerator = df_numerator,
    df_denominator = df_denominator,
    alpha = res$alpha,
    cv_score = switch(cv,
      res$loss,
      NULL
    ),
    scale = scale,
    sigma = sigma,
    centers = df_centers,
    model_matrices = dat,
    nfold = switch(cv,
      max(cv_ind_nu) + 1,
      NULL
    ),
    constrained = constrained,
    alpha_opt = switch(cv,
      res$alpha[, which.min(res$loss), drop = FALSE],
      NULL
    ),
    sigma_opt = switch(cv,
      sigma[which.min(res$loss)],
      NULL
    ),
    call = cl
  )

  class(out) <- "kmm"
  out
}
