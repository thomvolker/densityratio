#' Single permutation
#' @rdname permute
#' @param object Object of class \code{ulsif} or \code{kliep}
#' @param ... Additional arguments to pass through to specific permute functions.
#' @return permutation statistic for a single permutation of the data
#' @importFrom stats update
#' @importFrom stats predict

permute <- function(object, ...) {
  UseMethod("permute")
}


#' Single permutation statistic of \code{ulsif} object
#' @rdname permute
#' @param object \code{ulsif} object
#' @param stacked \code{matrix} with stacked numerator and denominator samples
#' @param nnu Scalar with numerator sample size
#' @param nde Scalar with denominator sample size
#' @return permutation statistic for a single permutation of the data
#' @method permute ulsif
#' @importFrom stats predict
#' @importFrom stats update

permute.ulsif <- function(object, stacked, nnu, nde) {
  ind <- sample(rep(c(TRUE, FALSE), times = c(nnu, nde)))
  r <- stats::update(
    object,
    df_numerator = stacked[ind, ],
    df_denominator = stacked[!ind, ],
    progressbar = FALSE
  )
  1/(2 * nnu) * sum(c(stats::predict(r, newdata = stacked[ind, ]))) - 1/nde * sum(c(stats::predict(r, newdata = stacked[!ind, ]))) + 1/2
}
