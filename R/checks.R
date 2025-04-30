#' @importFrom stats quantile
#' @importFrom stats median
#' @importFrom stats model.matrix
#' @importFrom stats sd

check.datatype <- function(data) {
  if (is.vector(data)) {
    data <- data.frame(x = data)
  } else {
    data <- as.data.frame(data)
  }

  if (sum(is.na(data)) > 0) {
    stop("Missing data can currently not be handled, please solve the missing data problem first.")
  }
  data
}

check.dataform <- function(nu, de, centers = NULL, nullcenters, newdata = NULL, scale) {
  numvars <- which(sapply(nu, is.numeric))
  if (!is.null(centers)) {
    alldat <- rbind(nu, de, centers)
    ind <- c(rep("nu", nrow(nu)), rep("de", nrow(de)), rep("ce", nrow(centers)))
  } else {
    alldat <- rbind(nu, de)
    ind <- c(rep("nu", nrow(nu)), rep("de", nrow(de)))
  }

  if (!is.null(scale)) {
    scale <- match.arg(scale, c("numerator", "denominator"))
  }

  scaledat <- if (identical(scale, "numerator")) {
    nu[, numvars, drop = FALSE]
  } else if (identical(scale, "denominator")) {
    de[, numvars, drop = FALSE]
  }

  if (!is.null(scale)) {
    if (!nullcenters) {
      warning("Note that you provided centers while also applying scaling to the variables. The centers are scaled accordingly.")
    }

    means <- colMeans(scaledat)
    sds <- sapply(scaledat, sd)

    alldat[, numvars] <- scale(alldat[, numvars], center = means, scale = sds) |> as.data.frame()

    if (any(sds == 0)) {
      warning("Some variables have zero variance in the data used for scaling. These variables are removed from both the numerator and denominator data.")
      remove_vars <- numvars[which(sds == 0)]
      alldat <- alldat[, -remove_vars, drop = FALSE]
    }
  }

  dat <- model.matrix(~ . - 1, alldat)

  if (!is.null(newdata)) {
    newdata <- check.datatype(newdata)
    if (!is.null(scale)) {
      newdata[, numvars] <- scale(newdata[, numvars], center = means, scale = sds) |> as.data.frame()
    }
    alldat <- rbind(alldat, newdata)
    newdat <- model.matrix(~ . - 1, alldat)
    newdat <- newdat[(nrow(dat) + 1):(nrow(dat) + nrow(newdata)), , drop = FALSE]
    dat <- newdat[, colnames(dat), drop = FALSE]
    return(dat)
  } else {
    dat <- list(
      nu = dat[ind == "nu", , drop = FALSE],
      de = dat[ind == "de", , drop = FALSE],
      ce = dat[ind == "ce", , drop = FALSE]
    )
    return(dat)
  }
}

check.variables <- function(nu, de, ce = NULL) {
  nu <- check.datatype(nu)
  de <- check.datatype(de)

  numvars_nu <- which(sapply(nu, is.numeric))
  numvars_de <- which(sapply(de, is.numeric))

  if (
    !identical(numvars_nu, numvars_de) |
      !identical(ncol(nu), ncol(de)) |
      !identical(colnames(nu), colnames(de))
  ) {
    stop("The numerator and denominator data must contain the same variables.")
  }
  if (!is.null(ce)) {
    ce <- check.datatype(ce)
    numvars_ce <- which(sapply(ce, is.numeric))
    if (
      !identical(numvars_nu, numvars_ce) |
        !identical(ncol(nu), ncol(ce)) |
        !identical(colnames(nu), colnames(ce))
    ) {
      stop("The data and centers must contain the same variables.")
    }
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
      sigma <- stats::quantile(dist_nu, p)
      sigma <- sqrt(sigma / 2)
    }
  }
  # if both sigma and sigma_quantile are not specified, specify the quantiles linearly, based on nsigma
  else {
    if (!is.numeric(nsigma) | length(nsigma) != 1) {
      stop("If 'sigma_quantile' and 'sigma' are not specified, 'nsigma' must be a scalar value.")
    } else {
      if (nsigma < 1) {
        stop("'nsigma' must be a positive scalar.")
      } else if (nsigma == 1) {
        sigma <- stats::median(dist_nu)
        sigma <- sqrt(sigma / 2)
      } else {
        p <- seq(0.05, 0.95, length.out = nsigma)
        sigma <- stats::quantile(dist_nu, p)
        sigma <- sqrt(sigma / 2)
      }
    }
  }
  u <- unique(sigma)
  if (length(u) < length(sigma)) {
    warning("There are duplicate values in 'sigma', only the unique values are used.\n")
  }
  if (any(u <= 0)) {
    stop("The values in 'sigma' must be larger than 0.")
  }
  u
}

check.sigma_quantile.lhss <- function(nsigma, sigma, sigma_quantile) {
  if (!is.null(sigma)) {
    if (!is.numeric(sigma) | !is.null(dim(sigma))) {
      stop("If 'sigma' is specified it must be a numeric scalar or vector.")
    } else {
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
      } else if (nsigma == 1) {
        p <- 0.5
      } else {
        p <- seq(0.05, 0.95, length.out = nsigma)
      }
    }
  }
  p
}

check.lambda <- function(nlambda, lambda) {
  if (!is.null(lambda)) {
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

check.centers <- function(dat, centers, ncenters) {
  if (!is.null(centers)) {
    centers <- check.datatype(centers)
  } else {
    if (!is.numeric(ncenters) | length(ncenters) != 1 | ncenters < 1) {
      stop("The 'ncenters' parameter must be a positive numeric scalar.")
    } else if (ncenters >= nrow(dat)) {
      centers <- dat
    } else {
      centers <- dat[sample(nrow(dat), ncenters), , drop = FALSE]
    }
  }
  centers
}

check.intercept <- function(intercept) {
  if (!(isTRUE(intercept) | isFALSE(intercept))) {
    stop("'intercept' must be either 'TRUE' or 'FALSE'")
  }
  intercept
}

check.symmetric <- function(nu, centers) {
  if (isTRUE(all.equal(nu, centers, check.attributes = FALSE))) {
    symmetric <- TRUE
  } else {
    symmetric <- FALSE
  }
  symmetric
}

check.parallel <- function(parallel, nthreads, iterator) {
  if (!(isTRUE(parallel) | isFALSE(parallel))) {
    stop("'parallel' must be either 'TRUE' or 'FALSE'")
  }
  if (parallel & length(iterator) == 1) {
    warning("No argument to parallelize over.\n")
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
      if (!is.numeric(nthreads) | !(length(nthreads) == 1)) {
        stop("If parallel processing is enabled, nthreads must be NULL or a positive integer.")
      }
      if (nthreads < 1) {
        nthreads <- 1
        warning("The minimum number of threads equals 1, 'nthreads' is set to 1.\n")
      }
    }
  }
  nthreads
}

check.epsilon <- function(epsilon) {
  if (!(is.null(epsilon))) {
    if (!is.numeric(epsilon) | !is.vector(epsilon) | any(epsilon < 0)) {
      stop("'eps' must be either NULL or a positive scalar or vector.")
    }
  } else {
    epsilon <- 10^{
      1:-5
    }
  }
  epsilon
}

check.maxit <- function(maxit) {
  if (maxit < 0 | !is.numeric(maxit) | length(maxit) != 1 | !is.finite(maxit)) {
    stop("'maxit' must be a positive scalar")
  }
  maxit
}

check.nfold <- function(cv, nfold, n) {
  if (!cv) {
    cv_ind <- rep(0, n)
  } else {
    if (!is.numeric(nfold) | length(nfold) > 1) {
      stop("'nfold' must be a positive scalar indicating the number of cross-validation folds.")
    } else if (nfold < 2 | nfold > n) {
      stop(paste0("'nfold' must be at least 2 and at most ", n, ", given the number of observations"))
    } else {
      cv_ind <- sample(rep_len(0:(nfold - 1), n))
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
      } else if (length(object$sigma) == 1) {
        sigma <- object$sigma
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

    if (sigma == "sigmaopt") {
      # for all lambda values considered, select optimal if it exists, NA otherwise
      sigmaind <- sapply(lambdaind, function(i) {
        ifelse(is.na(i), NA, which.min(object$cv_score[, i]))
      })
      lambdasigmaind <- matrix(c(lambdaind, sigmaind), ncol = 2)
    } else if (sigma == "all") {
      # for all lambda values considered, select all sigma values if they already
      # exist, NA otherwise
      lambdasigmaind <- matrix(
        c(
          rep(lambdaind, each = nrow(object$sigma)),
          rep(seq_len(nrow(object$sigma)), nlambda)
        ),
        ncol = 2
      )
      lambdasigmaind[is.na(lambdasigmaind[,1]), 2] <- NA
    }
  } else if (is.numeric(sigma) & is.vector(sigma)) {
    # check if any of the lambda-sigma combinations are already fitted, and set
    # indices to NA otherwise
    sigmaind <- unlist(
      lapply(lambdaind, function(i) {
        if (is.na(i)) {
          rep(NA, length(sigma))
        } else {
          match(sigma, object$sigma[, i])
        }
      })
    )
    lambdasigmaind <- matrix(
      c(
        rep(lambdaind, each = length(sigma)),
        sigmaind
      ),
      ncol = 2
    )
  } else {
    stop("'sigma' should be one of 'sigmaopt', 'all' or a numeric scalar or vector with values to use as sigma parameter")
  }
  lambdanew <- rep(lambda, each = nrow(lambdasigmaind)/nlambda)
  sigmanew <- sapply(seq_len(nrow(lambdasigmaind)), function(i) {

    if (any(is.na(lambdasigmaind[i, 1:2]))) {
      if (is.numeric(sigma)) {
        return(sigma[(i - 1) %% length(sigma) + 1])
      } else {
        return(NA)
      }
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

check.subspace.spectral.predict <- function(object, m) {
  if (is.character(m)) {
    m <- match.arg(m, c("opt", "all"))
    if (m == "opt") {
      m <- object$m_opt
    } else if (m == "all") {
      m <- object$m
    }
  } else if (is.numeric(m) & is.vector(m)) {
    if (max(m) > nrow(object$centers)) {
      stop("'m' cannot exceed the number of centers.")
    }
  } else {
    stop("'m' should be one of 'opt', 'all' or an integer vector with values to use as 'm' parameter")
  }
  m
}

check.subspace <- function(m, P) {
  if (is.null(m)) {
    m <- floor(sqrt(P))
  } else {
    if (!is.numeric(m)) {
      stop("The dimension of the subspace 'm' must be 'NULL', an integer value, or an integer vector.")
    } else {
      if (any(m %% 1 != 0)) stop("The dimension of the subspace 'm' must be 'NULL', an integer value or an integer vector.")
      if (m > P) stop("The dimension of the subspace 'm' must be smaller than the number of variables.")
    }
  }
  m
}

check.subspace.spectral <- function(m, cv_ind_de) {
  ncenters <- ifelse(all(cv_ind_de == 0),
    length(cv_ind_de),
    length(cv_ind_de) - max(table(cv_ind_de))
  )

  if (is.null(m)) {
    m <- seq(1, ncenters, length.out = 50)
    m <- unique(floor(m))
  } else {
    if (!is.numeric(m)) {
      stop("The number of spectral components 'J' must be 'NULL' or a numeric scalar or vector.")
    } else {
      if (any(m < 1) | any(m > ncenters)) stop(paste0("The number of spectral components 'm' must be at least 1 and at most ", ncenters, "."))
      m <- unique(floor(m))
    }
  }
  m
}

check.newdata <- function(object, newdata) {
  if (!is.null(newdata)) {
    check.variables(
      check.datatype(object$df_numerator),
      check.datatype(newdata)
    )
    newdata <- check.dataform(
      check.datatype(object$df_numerator),
      check.datatype(object$df_denominator),
      NULL,
      TRUE,
      newdata,
      object$scale
    )
  } else {
    newdata <- object$model_matrices$nu
  }
  newdata
}

check.var.names <- function(vars, data) {
  nm <- colnames(data)
  if (!all(vars %in% nm)) {
    stop(
      paste0(
        "The indicated variables are not names of columns in the data. The variables are: ",
        paste0(nm, collapse = ", ")
      )
    )
  }
}

check.object.type <- function(object) {
  models <- c(
    "ulsif",
    "kliep",
    "kmm",
    "lhss",
    "spectral",
    "naivedensityratio",
    "naivesubspacedensityratio"
  )

  if (!attr(object, "class") %in% models) {
    stop("Only densityratio objects can be plotted.")
  }
}

check.logscale <- function(ext, logscale, tol) {
  if (logscale) {
    # Convert values lower than tolerance to tol
    negdr <- ext$dr < tol
    ext$dr[negdr] <- tol
    if (any(negdr)) {
      warning(
        paste(
          "Negative estimated density ratios for", sum(negdr),
          "observation(s) converted to", tol,
          "before applying logarithmic transformation"
        ),
        call. = FALSE
      )
    }

    # Apply log transformation
    ext$dr <- log(ext$dr)

    # Set y axis intercept to 0
    ext$yintercept <- 0
  } else {
    # Set y axis intercept to 1
    ext$yintercept <- 1
  }
  return(ext)
}
