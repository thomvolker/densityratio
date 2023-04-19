
check.dataform <- function(nu, de) {
  if (! (is.numeric(nu) & is.numeric(de))) {
    stop("Currently only numeric data is supported.")
  }
  if (ncol(nu) != ncol(de) | !all(names(nu) == names(de))) {
    stop("nu and de must contain exactly the same set of variables.")
  }
}

check.sigma <- function(sigma) {
  if (!(is.null(sigma) | is.numeric(sigma)) | !is.null(dim(sigma))) {
    stop("sigma must be either NULL, a numeric scalar or a numeric vector")
  }
}


check.lambda <- function(lambda) {
  if (!(is.null(lambda) | is.numeric(lambda)) | !is.null(dim(lambda))) {
    stop("lambda must be either NULL, a numeric scalar or a numeric vector")
  }
}

check.centers <- function(nu, centers, ncenters) {

  if (!is.null(centers)) {

    centers <- as.matrix(centers)
    check.dataform(nu, centers)

    if (!is.null(ncenters)) {
      if (nrow(centers) > ncenters) {
        warning("You have specified more centers than you allowed by ncenters, the ncenters parameter is ignored.")
      }
      if (nrow(centers) < ncenters) {
        warning("ncenters > nrow(centers), the numerator samples are used as gaussian centers.")
      }
    }
  }
  if (is.null(centers) & (is.null(ncenters) | ncenters == nrow(nu))) {
    centers <- nu
  }
  if (is.null(centers) & !is.null(ncenters)) {
    if (!is.numeric(ncenters) | length(ncenters) > 1) {
      stop("If ncenters is specified, it must be a scalar value.")
    } else {
      centers <- nu[sample(nrow(nu), ncenters), , drop = FALSE]
    }
  }

  centers
}

check.symmetric <- function(centers, ncenters) {
  if (is.null(centers) & (is.null(ncenters) | ncenters >= nrow(centers))) {
    TRUE
  } else {
    FALSE
  }
}

check.parallel <- function(parallel, nthreads, sigma, lambda) {
  if (!is.logical(parallel)) {
    stop("'parallel' must be either 'TRUE' or 'FALSE'")
  }
  if (parallel & length(sigma) == 1 & length(lambda) == 1) {
    warning("Parallel processing only works for multiple 'sigma' and/or 'lambda' values")
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
  }
  else {
    if(!is.null(nthreads)) {
      if(!(is.numeric(nthreads) & length(nthreads) == 1)) {
        stop("If parallel processing is enabled, nthreads must be NULL or a positive integer.")
      }
      if(nthreads < 1) {
        nthreads <- -1
      }
    } else {
      nthreads <- 0
    }
  }
  nthreads
}
