#' Single permutation
#' @rdname permute
#' @param object Density ratio object
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

permute.ulsif <- function(object, stacked, nnu, nde, ...) {
  ind <- sample(rep(c(TRUE, FALSE), times = c(nnu, nde)))
  r <- stats::update(
    object = object,
    df_numerator = stacked[ind, ],
    df_denominator = stacked[!ind, ],
    progressbar = FALSE
  )

  pred_nu <- c(stats::predict(r, newdata = stacked[ind, , drop = FALSE]))
  pred_de <- c(stats::predict(r, newdata = stacked[!ind, , drop = FALSE]))

  1/(2 * nnu) * sum(pred_nu) - 1/nde * sum(pred_de) + 1/2
}

#' Single permutation statistic of \code{kliep} object
#' @rdname permute
#' @param object \code{kliep} object
#' @param stacked \code{matrix} with stacked numerator and denominator samples
#' @param nnu Scalar with numerator sample size
#' @param nde Scalar with denominator sample size
#' @param min_pred Minimum value of the density ratio
#' @return permutation statistic for a single permutation of the data
#' @method permute kliep
#' @importFrom stats predict
#' @importFrom stats update

permute.kliep <- function(object, stacked, nnu, nde, min_pred = sqrt(.Machine$double.eps), ...) {
  ind <- sample(rep(c(TRUE, FALSE), times = c(nnu, nde)))
  r <- stats::update(
    object = object,
    df_numerator = stacked[ind, ],
    df_denominator = stacked[!ind, ],
    progressbar = FALSE
  )

  pred_nu <- c(stats::predict(r, newdata = stacked[ind, , drop = FALSE]))

  mean(log(pmax(min_pred, pred_nu)))
}

#' Single permutation statistic of \code{kmm} object
#' @rdname permute
#' @param object \code{kmm} object
#' @param stacked \code{matrix} with stacked numerator and denominator samples
#' @param nnu Scalar with numerator sample size
#' @param nde Scalar with denominator sample size
#' @return permutation statistic for a single permutation of the data
#' @method permute kmm
#' @importFrom stats predict
#' @importFrom stats update

permute.kmm <- function(object, stacked, nnu, nde, ...) {
  ind <- sample(rep(c(TRUE, FALSE), times = c(nnu, nde)))
  r <- stats::update(
    object = object,
    df_numerator = stacked[ind, ],
    df_denominator = stacked[!ind, ],
    progressbar = FALSE
  )

  pred_nu <- c(stats::predict(r, newdata = stacked[ind, , drop = FALSE]))
  pred_de <- c(stats::predict(r, newdata = stacked[!ind, , drop = FALSE]))

  1/(2 * nnu) * sum(pred_nu) - 1/nde * sum(pred_de) + 1/2
}

#' Single permutation statistic of \code{lhss} object
#' @rdname permute
#' @param object \code{lhss} object
#' @param stacked \code{matrix} with stacked numerator and denominator samples
#' @param nnu Scalar with numerator sample size
#' @param nde Scalar with denominator sample size
#' @return permutation statistic for a single permutation of the data
#' @method permute lhss
#' @importFrom stats predict
#' @importFrom stats update

permute.lhss <- function(object, stacked, nnu, nde, ...) {
  ind <- sample(rep(c(TRUE, FALSE), times = c(nnu, nde)))
  r <- stats::update(
    object = object,
    df_numerator = stacked[ind, ],
    df_denominator = stacked[!ind, ],
    progressbar = FALSE
  )

  pred_nu <- c(stats::predict(r, newdata = stacked[ind, , drop = FALSE]))
  pred_de <- c(stats::predict(r, newdata = stacked[!ind, , drop = FALSE]))

  1/(2 * nnu) * sum(pred_nu) - 1/nde * sum(pred_de) + 1/2
}

#' Single permutation statistic of \code{spectral} object
#' @rdname permute
#' @param object \code{spectral} object
#' @param stacked \code{matrix} with stacked numerator and denominator samples
#' @param nnu Scalar with numerator sample size
#' @param nde Scalar with denominator sample size
#' @return permutation statistic for a single permutation of the data
#' @method permute spectral
#' @importFrom stats predict
#' @importFrom stats update

permute.spectral <- function(object, stacked, nnu, nde, ...) {
  ind <- sample(rep(c(TRUE, FALSE), times = c(nnu, nde)))
  r <- stats::update(
    object = object,
    df_numerator = stacked[ind, ],
    df_denominator = stacked[!ind, ],
    progressbar = FALSE
  )

  pred_nu <- c(stats::predict(r, newdata = stacked[ind, , drop = FALSE]))
  pred_de <- c(stats::predict(r, newdata = stacked[!ind, , drop = FALSE]))

  1/(2 * nnu) * sum(pred_nu) - 1/nde * sum(pred_de) + 1/2
}

#' Single permutation statistic of \code{naivedensityratio} object
#' @rdname permute
#' @param object \code{naivedensityratio} object
#' @param stacked \code{matrix} with stacked numerator and denominator samples
#' @param nnu Scalar with numerator sample size
#' @param nde Scalar with denominator sample size
#' @param min_pred Minimum value of the predicted density ratio
#' @param max_pred Maximum value of the predicted density ratio
#' @return permutation statistic for a single permutation of the data
#' @method permute naivedensityratio
#' @importFrom stats predict
#' @importFrom stats update

permute.naivedensityratio <- function(object, stacked, nnu, nde, min_pred, max_pred) {
  ind <- sample(rep(c(TRUE, FALSE), times = c(nnu, nde)))
  r <- stats::update(
    object = object,
    df_numerator = stacked[ind, ],
    df_denominator = stacked[!ind, ],
  )

  pred_nu <- c(stats::predict(r, newdata = stacked[ind, , drop = FALSE]))
  pred_de <- c(stats::predict(r, newdata = stacked[!ind, , drop = FALSE]))

  (mean(pmin(max_pred, pmax(min_pred, pred_nu))) - mean(pmin(max_pred, pmax(min_pred, pred_de))))^2
}


