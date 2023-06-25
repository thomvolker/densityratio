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
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @export


summary.ulsif <- function(object, test = TRUE, n.perm = 100, parallel = FALSE, cl = NULL, ...) {

  nu <- as.matrix(object$df_numerator)
  de <- as.matrix(object$df_denominator)
  stacked <- rbind(nu, de)
  nnu <- nrow(nu)
  nde <- nrow(de)

  if (parallel & is.null(cl)) {
    ncores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(ncores)
  }

  pred_nu <- c(stats::predict(object, newdata = nu))
  pred_de <- c(stats::predict(object, newdata = de))

  obsPE  <- 1/(2 * nnu) * sum(c(stats::predict(object, newdata = nu))) - 1/nde * sum(c(stats::predict(object, newdata = de))) + 1/2
  if (test) {
    distPE <- pbapply::pbreplicate(
      n.perm,
      permute(object, stacked = stacked, nnu = nnu, nde = nde),
      simplify = TRUE,
      cl = cl
    )
  }

  out <- list(
    alpha_opt  = object$alpha_opt,
    sigma_opt  = object$sigma_opt,
    lambda_opt = object$lambda_opt,
    centers    = object$centers,
    dr = data.frame(dr = c(pred_nu, pred_de),
                    group = factor(c(rep(c("nu", "de"), c(nnu, nde))))),
    PE = obsPE,
    refPE = switch(test, distPE, NULL),
    p_value = switch(test, mean(distPE > obsPE), NULL),
    call = object$call
  )
  class(out) <- "summary.ulsif"
  out
}


summary.ulsif <- function(object, test = TRUE, n.perm = 100, parallel = FALSE, cl = NULL, ...) {

  nu <- as.matrix(object$df_numerator)
  de <- as.matrix(object$df_denominator)
  stacked <- rbind(nu, de)
  nnu <- nrow(nu)
  nde <- nrow(de)

  if (parallel & is.null(cl)) {
    ncores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(ncores)
  }

  pred_nu <- c(stats::predict(object, newdata = nu))
  pred_de <- c(stats::predict(object, newdata = de))

  obsPE  <- 1/(2 * nnu) * sum(c(stats::predict(object, newdata = nu))) - 1/nde * sum(c(stats::predict(object, newdata = de))) + 1/2
  if (test) {
    distPE <- pbapply::pbreplicate(
      n.perm,
      permute(object, stacked = stacked, nnu = nnu, nde = nde),
      simplify = TRUE,
      cl = cl
    )
  }

  out <- list(
    alpha_opt  = object$alpha_opt,
    sigma_opt  = object$sigma_opt,
    lambda_opt = object$lambda_opt,
    centers    = object$centers,
    dr = data.frame(dr = c(pred_nu, pred_de),
                    group = factor(c(rep(c("nu", "de"), c(nnu, nde))))),
    PE = obsPE,
    refPE = switch(test, distPE, NULL),
    p_value = switch(test, mean(distPE > obsPE), NULL),
    call = object$call
  )
  class(out) <- "summary.ulsif"
  out
}
