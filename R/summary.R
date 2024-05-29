#' Extract summary from \code{ulsif} object, including two-sample significance
#' test for homogeneity of the numerator and denominator samples
#'
#' @rdname summary.ulsif
#' @param object Object of class \code{ulsif}
#' @param test logical indicating whether to statistically test for homogeneity
#' of the numerator and denominator samples.
#' @param n_perm Scalar indicating number of permutation samples
#' @param parallel \code{logical} indicating to run the permutation test in parallel
#' @param cl \code{NULL} or a cluster object created by \code{makeCluster}. If
#' \code{NULL} and \code{parallel = TRUE}, it uses the number of available cores
#' minus 1.
#' @param ... further arguments passed to or from other methods.
#' @return Summary of the fitted density ratio model
#' @method summary ulsif
#' @importFrom stats predict
#' @importFrom pbapply pbreplicate
#' @importFrom parallel detectCores makeCluster stopCluster
#' @export


summary.ulsif <- function(object, test = FALSE, n_perm = 100, parallel = FALSE, cl = NULL, ...) {

  nu <- check.datatype(object$df_numerator)
  de <- check.datatype(object$df_denominator)
  stacked <- rbind(nu, de)
  nnu <- nrow(nu)
  nde <- nrow(de)

  if (parallel & is.null(cl)) {
    ncores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl))
  }

  pred_nu <- c(stats::predict(object, newdata = object$df_numerator))
  pred_de <- c(stats::predict(object, newdata = object$df_denominator))

  obsPE  <- 1/(2 * nnu) * sum(pred_nu) - 1/nde * sum(pred_de) + 1/2
  if (test) {
    distPE <- pbapply::pbreplicate(
      n_perm,
      permute(object, stacked = stacked, nnu = nnu, nde = nde),
      simplify = TRUE,
      cl = cl
    )
    rec <- update(object, df_numerator = de, df_denominator = nu, progressbar = FALSE)
    recPE <- 1/(2 * nde) * sum(c(stats::predict(rec, newdata = de))) -
      1/nnu * sum(c(stats::predict(rec, newdata = nu))) + 1/2
  }

  out <- list(
    alpha_opt  = object$alpha_opt,
    sigma_opt  = object$sigma_opt,
    lambda_opt = object$lambda_opt,
    centers    = object$centers,
    dr = data.frame(dr = c(pred_nu, pred_de),
                    group = factor(rep(c("nu", "de"), c(nnu, nde)))),
    PE = obsPE,
    PErec = switch(test, recPE, NULL),
    refPE = switch(test, distPE, NULL),
    p_value = switch(test, min(1, 2 * mean(distPE > max(obsPE, recPE))), NULL),
    call = object$call
  )
  class(out) <- "summary.ulsif"
  out
}


#' Extract summary from \code{kliep} object, including two-sample significance
#' test for homogeneity of the numerator and denominator samples
#'
#' @rdname summary.kliep
#' @param object Object of class \code{kliep}
#' @param test logical indicating whether to statistically test for homogeneity
#' of the numerator and denominator samples.
#' @param n_perm Scalar indicating number of permutation samples
#' @param parallel \code{logical} indicating to run the permutation test in parallel
#' @param cl \code{NULL} or a cluster object created by \code{makeCluster}. If
#' \code{NULL} and \code{parallel = TRUE}, it uses the number of available cores
#' minus 1.
#' @param min_pred Scalar indicating the minimum value for the predicted density
#' ratio values (used in the divergence statistic) to avoid negative density
#' ratio values.
#' @param ... further arguments passed to or from other methods.
#' @return Summary of the fitted density ratio model
#' @method summary kliep
#' @importFrom stats predict
#' @importFrom pbapply pbreplicate
#' @importFrom parallel detectCores makeCluster stopCluster
#' @export

summary.kliep <- function(object, test = FALSE, n_perm = 100, parallel = FALSE, cl = NULL, min_pred = 1e-6, ...) {

  nu <- check.datatype(object$df_numerator)
  de <- check.datatype(object$df_denominator)
  stacked <- rbind(nu, de)
  nnu <- nrow(nu)
  nde <- nrow(de)

  if (parallel & is.null(cl)) {
    ncores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl))
  }

  pred_nu <- c(stats::predict(object, newdata = nu))
  pred_de <- c(stats::predict(object, newdata = de))

  obsUKL  <- mean(log(pmax(min_pred, c(pred_nu))))
  if (test) {
    distUKL <- pbapply::pbreplicate(
      n_perm,
      permute(object, stacked = stacked, nnu = nnu, nde = nde),
      simplify = TRUE,
      cl = cl
    )
    rec <- update(object, df_numerator = de, df_denominator = nu, progressbar = FALSE)
    recUKL <- mean(log(pmax(sqrt(.Machine$double.eps), c(stats::predict(rec, newdata = de)))))
  }

  out <- list(
    alpha_opt  = object$alpha_opt,
    sigma_opt  = object$sigma_opt,
    centers    = object$centers,
    dr = data.frame(dr = c(pred_nu, pred_de),
                    group = factor(rep(c("nu", "de"), c(nnu, nde)))),
    UKL = obsUKL,
    UKLrec = switch(test, recUKL, NULL),
    refUKL = switch(test, distUKL, NULL),
    p_value = switch(test, min(1, 2 * mean(distUKL > max(obsUKL, recUKL))), NULL),
    call = object$call
  )
  class(out) <- "summary.kliep"
  out
}

#' Extract summary from \code{lhss} object, including two-sample significance
#' test for homogeneity of the numerator and denominator samples
#'
#' @rdname summary.lhss
#' @param object Object of class \code{lhss}
#' @param test logical indicating whether to statistically test for homogeneity
#' of the numerator and denominator samples.
#' @param n_perm Scalar indicating number of permutation samples
#' @param parallel \code{logical} indicating to run the permutation test in parallel
#' @param cl \code{NULL} or a cluster object created by \code{makeCluster}. If
#' \code{NULL} and \code{parallel = TRUE}, it uses the number of available cores
#' minus 1.
#' @param ... further arguments passed to or from other methods.
#' @return Summary of the fitted density ratio model
#' @method summary lhss
#' @importFrom stats predict
#' @importFrom pbapply pbreplicate
#' @importFrom parallel detectCores makeCluster stopCluster
#' @export

summary.lhss <- function(object, test = FALSE, n_perm = 100, parallel = FALSE, cl = NULL, ...) {

  nu <- check.datatype(object$df_numerator)
  de <- check.datatype(object$df_denominator)
  stacked <- rbind(nu, de)
  nnu <- nrow(nu)
  nde <- nrow(de)

  if (parallel & is.null(cl)) {
    ncores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl))
  }

  pred_nu <- c(stats::predict(object, newdata = nu))
  pred_de <- c(stats::predict(object, newdata = de))

  obsPE  <- 1/(2 * nnu) * sum(pred_nu) - 1/nde * sum(pred_de) + 1/2
  if (test) {
    distPE <- pbapply::pbreplicate(
      n_perm,
      permute(object, stacked = stacked, nnu = nnu, nde = nde),
      simplify = TRUE,
      cl = cl
    )
    rec <- update(object, df_numerator = de, df_denominator = nu, progressbar = FALSE)
    recPE <- 1/(2 * nde) * sum(c(stats::predict(rec, newdata = de))) -
      1/nnu * sum(c(stats::predict(rec, newdata = nu))) + 1/2
  }

  out <- list(
    alpha_opt  = object$alpha_opt,
    sigma_opt  = object$sigma_opt,
    lambda_opt = object$lambda_opt,
    m = object$m,
    centers    = object$centers,
    dr = data.frame(dr = c(pred_nu, pred_de),
                    group = factor(rep(c("nu", "de"), c(nnu, nde)))),
    PE = obsPE,
    PErec = switch(test, recPE, NULL),
    refPE = switch(test, distPE, NULL),
    p_value = switch(test, min(1, 2 * mean(distPE > max(obsPE, recPE))), NULL),
    call = object$call
  )
  class(out) <- "summary.lhss"
  out
}

#' Extract summary from \code{spectral} object, including two-sample significance
#' test for homogeneity of the numerator and denominator samples
#' @rdname summary.spectral
#' @param object Object of class \code{spectral}
#' @param test logical indicating whether to statistically test for homogeneity
#' of the numerator and denominator samples.
#' @param n_perm Scalar indicating number of permutation samples
#' @param parallel \code{logical} indicating to run the permutation test in parallel
#' @param cl \code{NULL} or a cluster object created by \code{makeCluster}. If
#' \code{NULL} and \code{parallel = TRUE}, it uses the number of available cores
#' minus 1.
#' @param ... further arguments passed to or from other methods.
#' @return Summary of the fitted density ratio model
#' @method summary spectral
#' @importFrom stats predict
#' @importFrom pbapply pbreplicate
#' @importFrom parallel detectCores makeCluster stopCluster
#' @export


summary.spectral <- function(object, test = FALSE, n_perm = 100, parallel = FALSE, cl = NULL, ...) {

  nu <- check.datatype(object$df_numerator)
  de <- check.datatype(object$df_denominator)
  stacked <- rbind(nu, de)
  nnu <- nrow(nu)
  nde <- nrow(de)

  if (parallel & is.null(cl)) {
    ncores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl))
  }

  pred_nu <- c(stats::predict(object, newdata = nu, min_pred = sqrt(.Machine$double.eps)))
  pred_de <- c(stats::predict(object, newdata = de, min_pred = sqrt(.Machine$double.eps)))

  obsPE  <- 1/(2 * nnu) * sum(pred_nu) - 1/nde * sum(pred_de) + 1/2
  if (test) {
    distPE <- pbapply::pbreplicate(
      n_perm,
      permute(object, stacked = stacked, nnu = nnu, nde = nde, min_pred = sqrt(.Machine$double.eps)),
      simplify = TRUE,
      cl = cl
    )
    rec <- update(object, df_numerator = de, df_denominator = nu, progressbar = FALSE)
    recPE <- 1/(2 * nde) * sum(c(stats::predict(rec, newdata = de))) -
      1/nnu * sum(c(stats::predict(rec, newdata = nu))) + 1/2
  }

  out <- list(
    alpha_opt  = object$alpha_opt,
    sigma_opt  = object$sigma_opt,
    J_opt = object$J_opt,
    centers    = object$centers,
    dr = data.frame(dr = c(pred_nu, pred_de),
                    group = factor(rep(c("nu", "de"), c(nnu, nde)))),
    PE = obsPE,
    PErec = switch(test, recPE, NULL),
    refPE = switch(test, distPE, NULL),
    p_value = switch(test, min(1, 2 * mean(distPE > max(obsPE, recPE))), NULL),
    call = object$call
  )
  class(out) <- "summary.spectral"
  out
}

#' Extract summary from \code{naivedensityraito} object, including two-sample
#' significance test for homogeneity of the numerator and denominator samples
#'
#' @rdname summary.naivedensityratio
#' @param object Object of class \code{naivedensityratio}
#' @param test logical indicating whether to statistically test for homogeneity
#' of the numerator and denominator samples.
#' @param n_perm Scalar indicating number of permutation samples
#' @param parallel \code{logical} indicating to run the permutation test in parallel
#' @param cl \code{NULL} or a cluster object created by \code{makeCluster}. If
#' \code{NULL} and \code{parallel = TRUE}, it uses the number of available cores
#' minus 1.
#' @param ... further arguments passed to or from other methods.
#' @return Summary of the fitted density ratio model
#' @method summary naivedensityratio
#' @importFrom stats predict
#' @importFrom pbapply pbreplicate
#' @importFrom parallel detectCores makeCluster stopCluster
#' @export

summary.naivedensityratio <- function(object, test = FALSE, n_perm = 100, parallel = FALSE, cl = NULL, ...) {

  nu <- check.datatype(object$df_numerator)
  de <- check.datatype(object$df_denominator)

  stacked <- rbind(nu, de)

  nnu <- nrow(nu)
  nde <- nrow(de)

  if (parallel & is.null(cl)) {
    ncores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl))
  }

  pred_nu <- c(stats::predict(object, newdata = nu))
  pred_de <- c(stats::predict(object, newdata = de))

  min_pred <- log(sqrt(.Machine$double.eps))
  max_pred <- -min_pred
  SALDRD   <- (mean(pmin(max_pred, pmax(min_pred, pred_nu))) - mean(pmin(max_pred, pmax(min_pred, pred_de))))^2

  if (test) {
    distSALDRD <- pbapply::pbreplicate(
      n_perm,
      permute(object, stacked = stacked, nnu = nnu, nde = nde, min_pred = min_pred, max_pred = max_pred),
      simplify = TRUE,
      cl = cl
    )
  }

  out <- list(
    n = c(nnu = nnu, nde = nde),
    nvars = ncol(nu),
    dr = data.frame(dr = c(pred_nu, pred_de),
                    group = factor(rep(c("nu", "de"), c(nnu, nde)))),
    SALDRD = SALDRD,
    refSALDRD = switch(test, distSALDRD, NULL),
    p_value = switch(test, mean(distSALDRD > SALDRD), NULL),
    call = object$call
  )
  class(out) <- "summary.naivedensityratio"
  out
}

#' Extract summary from \code{naivesubspacedensityraito} object, including two-sample
#' significance test for homogeneity of the numerator and denominator samples
#'
#' @rdname summary.naivesubspacedensityratio
#' @param object Object of class \code{naivesubspacedensityratio}
#' @param test logical indicating whether to statistically test for homogeneity
#' of the numerator and denominator samples.
#' @param n_perm Scalar indicating number of permutation samples
#' @param parallel \code{logical} indicating to run the permutation test in parallel
#' @param cl \code{NULL} or a cluster object created by \code{makeCluster}. If
#' \code{NULL} and \code{parallel = TRUE}, it uses the number of available cores
#' minus 1.
#' @param ... further arguments passed to or from other methods.
#' @return Summary of the fitted density ratio model
#' @method summary naivesubspacedensityratio
#' @importFrom stats predict
#' @importFrom pbapply pbreplicate
#' @importFrom parallel detectCores makeCluster stopCluster
#' @export

summary.naivesubspacedensityratio <- function(object, test = FALSE, n_perm = 100, parallel = FALSE, cl = NULL, ...) {

  nu <- check.datatype(object$df_numerator)
  de <- check.datatype(object$df_denominator)

  stacked <- rbind(nu, de)

  nnu <- nrow(nu)
  nde <- nrow(de)

  if (parallel & is.null(cl)) {
    ncores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl))
  }

  pred_nu <- c(stats::predict(object, newdata = nu))
  pred_de <- c(stats::predict(object, newdata = de))

  min_pred <- log(sqrt(.Machine$double.eps))
  max_pred <- -min_pred
  SALDRD   <- (mean(pmin(max_pred, pmax(min_pred, pred_nu))) - mean(pmin(max_pred, pmax(min_pred, pred_de))))^2

  if (test) {
    distSALDRD <- pbapply::pbreplicate(
      n_perm,
      permute(object, stacked = stacked, nnu = nnu, nde = nde, min_pred = min_pred, max_pred = max_pred),
      simplify = TRUE,
      cl = cl
    )
  }

  out <- list(
    n = c(nnu = nnu, nde = nde),
    nvars = ncol(nu),
    subspace_dim = object$m,
    dr = data.frame(dr = c(pred_nu, pred_de),
                    group = factor(rep(c("nu", "de"), c(nnu, nde)))),
    SALDRD = SALDRD,
    refSALDRD = switch(test, distSALDRD, NULL),
    p_value = switch(test, mean(distSALDRD > SALDRD), NULL),
    call = object$call
  )
  class(out) <- "summary.naivesubspacedensityratio"
  out
}
