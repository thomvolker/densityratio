#' @importFrom stats quantile
#' @importFrom stats median
check.dataform <- function(nu, de) {
  # Only accept numeric input matrices (vectors/data.frames/tibbles
  # are converted to matrices by default)
  if (! (is.numeric(nu) & is.numeric(de))) {
    stop("Currently only numeric data is supported.")
  }
  # check whether the data sets have the same set of variables
  if (ncol(nu) != ncol(de) | !all.equal(colnames(nu), colnames(de))) {
    stop("nu and de must contain exactly the same set of variables.")
  }
  #
  if ((sum(is.na(nu)) + sum(is.na(de))) > 0) {
    stop("Your data has missing values, which cannot yet be handled.")
  }
}

check.sigma <- function(nsigma, sigma_quantile, sigma, dist_nu) {
  # disregard zero-distance observations
  dist_nu <- dist_nu[dist_nu > 0]
  # if sigma is specified, ignore 'nsigma' and 'sigma_quantile' and leave sigma as is
  if (!is.null(sigma)) {
    if (!is.numeric(sigma) | !is.null(dim(sigma))) {
      stop("If 'sigma' is specified it must be a numeric scalar or vector.")
    }
  }
  # if sigma is not specified, but sigma_quantile is, check the quantiles to be between 0 and 1
  # and ultimately set sigma to the specified quantiles of dist_nu
  else if (!is.null(sigma_quantile)) {
    if (!is.numeric(sigma_quantile) | !is.null(dim(sigma_quantile))) {
      stop("If 'sigma_quantile' is specified it must be a numeric scalar or vector.")
    } else if (min(sigma_quantile) <= 0 | max(sigma_quantile) >= 1) {
      stop("If 'sigma_quantile' is specified, the values must be larger than 0 and smaller than 1.")
    } else {
      p <- sigma_quantile
      sigma <- quantile(dist_nu, p) |> sqrt()
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
      else if (nsigma == 1) {
        sigma <- stats::median(dist_nu)
      } else {
        p <- seq(0.05, 0.95, length.out = nsigma)
        sigma <- quantile(dist_nu, p) |> sqrt()
      }
    }
  }
  u <- unique(sigma)
  if (length(u) < length(sigma)) {
    warning("There are duplicate values in 'sigma', only the unique values are used.\n")
  }
  u
}

check.sigma_quantile.lhss <- function(nsigma, sigma, sigma_quantile) {
  if (!is.null(sigma)) {
    if (!is.numeric(sigma) | !is.null(dim(sigma))) {
      stop("If 'sigma' is specified it must be a numeric scalar or vector.")
    } else{
      p <- NULL
    }
  } else if (!is.null(sigma_quantile)) {
    if (!is.numeric(sigma_quantile) | !is.null(dim(sigma_quantile))) {
      stop("If 'sigma_quantile' is specified it must be a numeric scalar or vector.")
    } else if (min(sigma_quantile) <= 0 | max(sigma_quantile) >= 1) {
      stop("If 'sigma_quantile' is specified, the values must be larger than 0 and smaller than 1.")
    } else {
      p <- sigma_quantile
    }
  } else {
    if (!is.numeric(nsigma) | length(nsigma) != 1) {
      stop("'nsigma' must be a scalar value.")
    } else {
      if (nsigma < 1) {
        stop("'nsigma' must be a positive scalar.")
      }
      else if (nsigma == 1) {
        p <- 0.5
      } else {
        p <- seq(0.05, 0.95, length.out = nsigma)
      }
    }
  }
  p
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
      warning(paste0("The 'ncenters' parameter exceeds the number of numerator records. Since the centers are selected from the numerator samples, 'ncenters' is set to ", nrow(nu), ".\n"))
    } else {
      centers <- nu[sample(nrow(nu), ncenters), , drop = FALSE]
    }
  }
  centers
}

check.intercept <- function(intercept) {
  if (!is.logical(intercept)) {
    stop("'intercept' must be either 'TRUE' or 'FALSE'")
  }
  intercept
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
    warning("Parallel processing only works for multiple 'lambda' values.\n")
    parallel <- FALSE
  }
  if (is.numeric(nthreads)) {
    if (parallel & nthreads < 2) {
      warning("Parallel processing only works for 2 or more threads; parallel processing is disabled.\n")
      parallel <- FALSE
    }
  }

  parallel

}

check.threads <- function(parallel, nthreads) {

  if (!parallel) {
    if (!is.null(nthreads)) {
      warning("Specifying 'nthreads' has no use without setting 'parallel = TRUE'.\n")
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
        warning("You specified a negative number of threads, 'nthreads' is set to 1.\n")
      }
    }
  }
  nthreads
}

check.epsilon <- function(epsilon) {
  if (!(is.null(epsilon))) {
    if (!is.numeric(epsilon) | !is.null(dim(epsilon))) {
      stop("'eps' must be either NULL, a numeric scalar or a numeric vector.")
    }
  } else {
    epsilon <- 10^{1:-5}
  }
  epsilon
}

check.maxit <- function(maxit) {
  if(maxit < 0) {
    stop("'maxit' must be a positive scalar")
  }
  maxit
}

check.nfold <- function(cv, nfold, sigma, n) {

  if (!cv) {
    cv_ind <- rep(0, n)
  } else {
    if (!is.numeric(nfold) | length(nfold) > 1) {
      stop("'nfold' must be a positive scalar indicating the number of cross-validation folds.")
    } else if (nfold < 2 | nfold > n) {
      stop(paste0("'nfold' must be at least 2 and at most ", n, ", given the number of observations"))
    } else {
      if (length(sigma) == 1) {
        warning("Only a single 'sigma' value is specified, while cross-validation serves to optimize 'sigma'.\n")
      }
      cv_ind <- sample(rep_len(0:(nfold-1), n))
    }
  }
  cv_ind
}

check.sigma.predict <- function(object, sigma) {
  if (is.character(sigma)) {
    sigma <- match.arg(sigma, c("sigmaopt", "all"))
    if (sigma == "sigmaopt") {
      if (!is.null(object$sigma_opt)) {
        sigma <- object$sigma_opt
      } else {
        warning("No optimal 'sigma' is defined, all 'sigma' values are used instead.")
        sigma <- object$sigma
      }
    } else if (sigma == "all") {
      sigma <- object$sigma
    }
  } else if (is.numeric(sigma) & is.vector(sigma)) {
    sigma <- sigma
  } else {
    stop("'sigma' should be one of 'sigmaopt', 'all' or a numeric scalar or vector with values to use as sigma parameter")
  }
  sigma
}

check.lambdasigma.predict <- function(object, sigma, lambda, lambdaind) {
  nlambda <- length(lambda)

  if (is.character(sigma)) {
    sigma <- match.arg(sigma, c("sigmaopt", "all"))
    if (sigma == "sigmaopt") { # Extract optimal sigma based on cv score for every lambda
      sigmaind <- sapply(lambdaind, \(i) which.min(object$cv_score[,i]))
      lambdasigmaind <- matrix(c(lambdaind, sigmaind), ncol = 2)
    } else if (sigma == "all") { # Extract all sigma values for every lambda
      sigmaind <- seq_len(nrow(object$sigma))
      lambdasigmaind <- matrix(c(rep(lambdaind, each = length(sigmaind)),
                                 rep(sigmaind, nlambda)),
                               ncol = 2)
      lambdasigmaind[is.na(lambdasigmaind[,1]), 2] <- NA
    }
  } else if (is.numeric(sigma) & is.vector(sigma)) {
    sigmaind <- lapply(lambdaind, \(i) match(sigma, object$sigma[,i]))
    lambdasigmaind <- matrix(c(rep(lambdaind, each = length(unlist(sigmaind))/length(lambdaind)),
                               unlist(sigmaind)),
                             ncol = 2)
  } else {
    stop("'sigma' should be one of 'sigmaopt', 'all' or a numeric scalar or vector with values to use as sigma parameter")
  }
  lambdanew <- rep(lambda, each = nrow(lambdasigmaind)/nlambda)
  sigmanew  <- sapply(1:nrow(lambdasigmaind), \(i) {
    if (is.numeric(sigma) & is.na(lambdasigmaind[i,2])) {
      return(sigma[(i-1) %% length(sigma) + 1])
    } else {
      return(object$sigma[lambdasigmaind[i, 2], lambdasigmaind[i, 1]])
    }
  })
  lambdasigma <- cbind(lambdasigmaind, lambdanew, sigmanew)
  colnames(lambdasigma) <- c("lambdaind", "sigmaind", "lambda", "sigma")
  lambdasigma
}

check.lambda.predict <- function(object, lambda) {
  if (is.character(lambda)) {
    lambda <- match.arg(lambda, c("lambdaopt", "all"))
    if (lambda == "lambdaopt") {
      lambda <- object$lambda_opt
    } else if (lambda == "all") {
      lambda <- object$lambda
    }
  } else if (is.numeric(lambda) & is.vector(lambda)) {
    lambda <- lambda
  } else {
    stop("'lambda' should be one of 'lambdaopt', 'all' or a numeric scalar or vector with values to use as lambda parameter")
  }
  lambda
}

check.J.predict <- function(object, J) {
  if (is.character(J)) {
    J <- match.arg(J, c("Jopt", "all"))
    if (J == "Jopt") {
      J <- object$J_opt
    } else if (J == "all") {
      J <- object$J
    }
  } else if (is.numeric(J) & is.vector(J)) {
    if (max(J) > nrow(object$centers)) {
      stop("'J' cannot exceed the number of centers.")
    }
  } else {
    stop("'J' should be one of 'Jopt', 'all' or an integer vector with values to use as J parameter")
  }
  J
}

check.subspace <- function(m, P) {
  if(is.null(m)) {
    m <- floor(sqrt(P))
  } else {
    if(!is.numeric(m)) {
      stop("The dimension of the subspace 'm' must be 'NULL', an integer value, or an integer vector.")
    } else {
      if (any(m %% 1 != 0)) stop("The dimension of the subspace 'm' must be 'NULL', an integer value or an integer vector.")
      if (m > P) stop("The dimension of the subspace 'm' must be smaller than the number of variables.")
    }
  }
  m
}

check.subspace.spectral <- function(J, cv_ind_de) {

  ncenters <- ifelse(all(cv_ind_de == 0),
                     length(cv_ind_de),
                     length(cv_ind_de) - max(table(cv_ind_de)))

  if (is.null(J)) {
    J <- seq(1, ncenters, length.out = 50)
    J <- unique(floor(J))
  } else {
    if (!is.numeric(J)) {
      stop("The number of spectral components 'J' must be 'NULL' or a numeric scalar or vector.")
    } else {
      if (any(J < 1) | any(J > ncenters)) stop(paste0("The number of spectral components 'J' must be at least 1 and at most ", ncenters, "."))
      J <- unique(floor(J))
    }
  }
  J
}

check.newdata <- function(object, newdata) {
  if (!is.null(newdata)) {
    newdata <- as.matrix(newdata)
    check.dataform(as.matrix(object$df_numerator), newdata)
  }
  else {
    newdata <- as.matrix(object$df_numerator)
  }
  newdata
}
