
check.dataform <- function(nu, de) {
  if (! (is.numeric(nu) & is.numeric(de))) {
    stop("Currently only numeric data is supported.")
  }
  if (ncol(nu) != ncol(de) | !all.equal(colnames(nu), colnames(de))) {
    stop("nu and de must contain exactly the same set of variables.")
  }
}

check.sigma <- function(nsigma, sigma_quantile, sigma, dist_nu) {

  # if sigma is specified, ignore 'nsigma' and 'sigma_quantile' and leave sigma as is
  if (!is.null(sigma)) {
    if (!is.numeric(sigma) | !is.null(dim(sigma))) {
      stop("If 'sigma' is specified it must be a numeric scalar or vector.")
    }
  }
  # if sigma is not specified, but sigma_quantile is, check the quantiles to be between 0 and 1
  # and ultimately set sigma to the specified quantiles of dist_nu
  else if (!is.null(sigma_quantile)) {
    if (min(sigma_quantile) <= 0 | max(sigma_quantile) >= 1) {
      stop("If 'sigma_quantile' is specified, the values must be larger than 0 and smaller than 1.")
    } else {
      p <- sigma_quantile
      sigma <- quantile(dist_nu, p)
    }
  }
  # if both sigma and sigma_quantile are not specified, specify the quantiles linearly, based on nsigma
  else {
    if (!is.numeric(nsigma) | length(nsigma) != 1) {
      stop("If 'sigma_quantile' and 'sigma' are not specified, 'nsigma' must be a scalar value.")
    } else {
      if (nsigma < 1) {
        stop("'nsigma' must be a positive scalar.")
      }
      p <- seq(1, nsigma) / (nsigma + 1)
      sigma <- quantile(dist_nu, p)
    }
  }
  sigma
}

check.lambda <- function(nlambda, lambda) {
  if (!(is.null(lambda))) {
    if (!is.numeric(lambda) | !is.null(dim(lambda))) {
      stop("'lambda' must be either NULL, a numeric scalar or a numeric vector.")
    }
  } else {
    if (!is.numeric(nlambda) | length(nlambda) != 1) {
      stop("'nlambda' must be a scalar")
    }
    lambda <- 10^seq(3, -3, length.out = nlambda)
  }
  lambda
}

check.centers <- function(nu, centers, ncenters) {

  if (!is.null(centers)) {
    centers <- as.matrix(centers)
    check.dataform(nu, centers)
  } else {
    if (!is.numeric(ncenters) | length(ncenters) != 1 | ncenters < 1) {
      stop("The 'ncenters' parameter must be a positive numeric scalar.")
    } else if (ncenters == nrow(nu)) {
      centers <- nu
    } else if (ncenters > nrow(nu)) {
      centers <- nu
      warning(paste0("The 'ncenters' parameter exceeds the number of numerator records. Since the centers are selected from the numerator samples, 'ncenters' is set to ", nrow(nu), "."))
    } else {
      centers <- nu[sample(nrow(nu), ncenters), , drop = FALSE]
    }
  }
  centers
}

check.symmetric <- function(nu, centers) {
  if(isTRUE(all.equal(nu, centers))) {
    symmetric <- TRUE
  } else {
    symmetric <- FALSE
  }
  symmetric
}

check.parallel <- function(parallel, nthreads, sigma, lambda) {
  if (!is.logical(parallel)) {
    stop("'parallel' must be either 'TRUE' or 'FALSE'")
  }
  if (parallel & length(lambda) == 1) {
    warning("Parallel processing only works for multiple 'lambda' values")
    parallel <- FALSE
  }
  if (is.numeric(nthreads)) {
    if (parallel & nthreads < 2) {
      warning("Parallel processing only works for 2 or more threads; parallel processing is disabled.")
      parallel <- FALSE
    }
  }

  parallel

}

check.threads <- function(parallel, nthreads) {

  if (!parallel) {
    if (!is.null(nthreads)) {
      warning("Specifying 'nthreads' has no use without setting 'parallel = TRUE'")
    }
    nthreads <- 0
  } else {
    if (is.null(nthreads)) {
      nthreads <- 0
    } else {
      if(!is.numeric(nthreads) | !(length(nthreads) == 1)) {
        stop("If parallel processing is enabled, nthreads must be NULL or a positive integer.")
      }
      if(nthreads < 1) {
        nthreads <- 1
        warning("You specified a negative number of threads, 'nthreads' is set to 1.")
      }
    }
  }
  nthreads
}
