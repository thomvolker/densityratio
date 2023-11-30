#' Extract summary from \code{ulsif} object, including two-sample significance
#' test for homogeneity of the numerator and denominator samples
#'
#' @rdname summary
#' @param object Object of class \code{ulsif} or \code{kliep}
#' @param test logical indicating whether to statistically test for homogeneity
#' of the numerator and denominator samples.
#' @param n.perm Scalar indicating number of permutation samples
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


summary.ulsif <- function(object, test = FALSE, n.perm = 100, parallel = FALSE, cl = NULL, ...) {

  nu <- as.matrix(object$df_numerator)
  de <- as.matrix(object$df_denominator)
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
      n.perm,
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
                    group = factor(c(rep(c("nu", "de"), c(nnu, nde))))),
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
#' @rdname summary
#' @return Summary of the fitted density ratio model
#' @method summary kliep
#' @importFrom stats predict
#' @importFrom pbapply pbreplicate
#' @importFrom parallel detectCores makeCluster stopCluster
#' @export

summary.kliep <- function(object, test = FALSE, n.perm = 100, parallel = FALSE, cl = NULL, ...) {

  nu <- as.matrix(object$df_numerator)
  de <- as.matrix(object$df_denominator)
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

  obsUKL  <- mean(log(pmax(sqrt(.Machine$double.eps), c(pred_nu))))
  if (test) {
    distUKL <- pbapply::pbreplicate(
      n.perm,
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

#' Extract summary from \code{naivedensityratio} object, including two-sample
#' significance test for homogeneity of the numerator and denominator samples
#'
#' @rdname summary
#' @return Summary of the fitted density ratio model
#' @method summary naivedensityratio
#' @importFrom stats predict
#' @importFrom pbapply pbreplicate
#' @importFrom parallel detectCores makeCluster stopCluster
#' @export

summary.naivedensityratio <- function(object, test = FALSE, n.perm = 100, parallel = FALSE, cl = NULL, ...) {

  nu <- as.matrix(object$df_numerator)
  de <- as.matrix(object$df_denominator)

  stacked <- rbind(nu, de)

  n_nu <- nrow(nu)
  n_de <- nrow(de)

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
      n.perm,
      permute(object, stacked = stacked, nnu = n_nu, nde = n_de, min_pred = min_pred, max_pred = max_pred),
      simplify = TRUE,
      cl = cl
    )
  }

  out <- list(
    n = c(n_nu = n_nu, n_de = n_de),
    nvars = ncol(nu),
    dr = data.frame(dr = c(pred_nu, pred_de),
                    group = factor(rep(c("nu", "de"), c(n_nu, n_de)))),
    SALDRD = SALDRD,
    refSALDRD = switch(test, distSALDRD, NULL),
    p_value = switch(test, mean(distSALDRD > SALDRD), NULL),
    call = object$call
  )
  class(out) <- "summary.naivedensityratio"
  out
}

#' Extract summary from \code{naivesubspacedensityratio} object, including
#' two-sample significance test for homogeneity of the numerator and
#' denominator samples
#'
#' @rdname summary
#' @return Summary of the fitted density ratio model
#' @method summary naivesubspacedensityratio
#' @importFrom stats predict
#' @importFrom pbapply pbreplicate
#' @importFrom parallel detectCores makeCluster stopCluster
#' @export

summary.naivesubspacedensityratio <- function(object, test = FALSE, n.perm = 100, parallel = FALSE, cl = NULL, ...) {

  nu <- as.matrix(object$df_numerator)
  de <- as.matrix(object$df_denominator)

  stacked <- rbind(nu, de)

  n_nu <- nrow(nu)
  n_de <- nrow(de)

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
      n.perm,
      permute(object, stacked = stacked, nnu = n_nu, nde = n_de, min_pred = min_pred, max_pred = max_pred),
      simplify = TRUE,
      cl = cl
    )
  }

  out <- list(
    n = c(n_nu = n_nu, n_de = n_de),
    nvars = ncol(nu),
    subspace_dim = object$m,
    dr = data.frame(dr = c(pred_nu, pred_de),
                    group = factor(rep(c("nu", "de"), c(n_nu, n_de)))),
    SALDRD = SALDRD,
    refSALDRD = switch(test, distSALDRD, NULL),
    p_value = switch(test, mean(distSALDRD > SALDRD), NULL),
    call = object$call
  )
  class(out) <- "summary.naivesubspacedensityratio"
  out
}
