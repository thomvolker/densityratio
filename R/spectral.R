#' Spectral series based density ratio estimation
#'
#' @param df_numerator \code{data.frame} with exclusively numeric variables with
#' the numerator samples
#' @param df_denominator \code{data.frame} with exclusively numeric variables
#' with the denominator samples (must have the same variables as
#' \code{df_denominator})
#' @param J Integer vector indicating the number of eigenvectors to use in the
#' spectral series expansion. Defaults to 50 evenly spaced values between 1 and
#' the number of denominator samples (or the largest number of samples that can
#' be used as centers in the cross-validation scheme).
#' @param scale \code{"numerator"}, \code{"denominator"}, or \code{FALSE},
#' indicating whether to standardize each numeric variable according to the
#' numerator means and standard deviations, the denominator means and standard
#' deviations, or apply no standardization at all.
#' @param nsigma Integer indicating the number of sigma values (bandwidth
#' parameter of the Gaussian kernel gram matrix) to use in cross-validation.
#' @param sigma_quantile \code{NULL} or numeric vector with probabilities to
#' calculate the quantiles of the distance matrix to obtain sigma values. If
#' \code{NULL}, \code{nsigma} values between \code{0.05} and \code{0.95} are
#' used.
#' @param sigma \code{NULL} or a scalar value to determine the bandwidth of the
#' Gaussian kernel gram matrix. If \code{NULL}, \code{nsigma} values between
#' \code{0.05} and \code{0.95} are used.
#' @param ncenters integer If smaller than the number of denominator observations,
#' an approximation to the eigenvector expansion based on only ncenters samples
#' is performed, instead of the full expansion. This can be useful for large
#' datasets.
#' @param cv logical indicating whether to use cross-validation to determine the
#' optimal sigma value and the optimal number of eigenvectors.
#' @param nfold Integer indicating the number of folds to use in the
#' cross-validation scheme. If \code{cv} is \code{FALSE}, this parameter is
#' ignored.
#' @param parallel logical indicating whether to use parallel processing in the
#' cross-validation scheme.
#' @param nthreads \code{NULL} or integer indicating the number of threads to
#' use for parallel processing. If parallel processing is enabled, it defaults
#' to the number of available threads minus one.
#' @param progressbar Logical indicating whether or not to display a progressbar.
#' @export
#'
#' @references Izbicki, R., Lee, A. &amp; Schafer, C. (2014). High-Dimensional
#' Density Ratio Estimation with Extensions to Approximate Likelihood Computation.
#' <i>Proceedings of Machine Learning Research</i> 33:420-429. Available from
#' <https://proceedings.mlr.press/v33/izbicki14.html>.


#' @return \code{spectral}-object, containing all information to calculate the
#' density ratio using optimal sigma and optimal spectral series expansion.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(100) |> matrix(100)
#' y <- rnorm(200, 1, 2) |> matrix(200)
#' spectral(x, y)
#' spectral(x, y, sigma = 2)

spectral <- function(df_numerator, df_denominator, J = NULL, scale = "numerator",
                     nsigma = 10, sigma_quantile = NULL, sigma = NULL,
                     ncenters = nrow(df_denominator), cv = TRUE, nfold = 10,
                     parallel = FALSE, nthreads = NULL, progressbar = TRUE) {

  cl  <- match.call()
  nu  <- check.datatype(df_numerator)
  de  <- check.datatype(df_denominator)

  check.variables(nu, de)
  nnu <- nrow(nu)
  nde <- nrow(de)

  ind <- if (ncenters < nde) sample(nde, ncenters) else seq_len(nde)
  df_centers <- de[ind, , drop = FALSE]
  dat <- check.dataform(nu, de, df_centers, TRUE, NULL, scale)

  #parallel  <- check.parallel(parallel, nthreads, sigma, lambda)
  nthreads  <- check.threads(parallel, nthreads)
  cv_ind_nu <- check.nfold(cv, nfold, sigma, nnu)
  cv_ind_de <- check.nfold(cv, nfold, sigma, ncenters)
  J         <- check.subspace.spectral(J, cv_ind_de)

  dist_nu <- distance(dat$nu, dat$ce, FALSE)
  dist_de <- distance(dat$ce, dat$ce, FALSE)

  sigma  <- check.sigma(nsigma, sigma_quantile, sigma, dist_nu)

  res <- spectral_dre(dist_nu, dist_de, J, sigma, cv_ind_nu, cv_ind_de,
                      parallel, nthreads, progressbar)

  # Order eigenvalues and eigenvectors from largest to smallest
  # note that armadillo orders them from smallest to largest
  res$Evals <- res$Evals[ncenters:(ncenters - max(J) + 1), , drop = FALSE]
  res$Evecs <- res$Evecs[, ncenters:(ncenters - max(J) + 1), , drop = FALSE]
  res$betatilde <- res$betatilde[max(J):1, , drop = FALSE]
  dimnames(res$Evecs) <- list(NULL, NULL, paste0("sigma", 1:length(sigma)))
  colnames(res$Evals) <- paste0("sigma", 1:length(sigma))

  colnames(res$betatilde) <- paste0("sigma", 1:length(sigma))

  opt_loss <- arrayInd(which.min(res$loss), dim(res$loss))

  out <- list(
    df_numerator = df_numerator,
    df_denominator = df_denominator,
    alpha = res$betatilde,
    Evecs = res$Evecs,
    Evals = res$Evals,
    cv_score = res$loss,
    sigma = sigma,
    J = J,
    centers = df_centers,
    model_matrices = dat,
    alpha_opt = res$betatilde[seq_len(J[opt_loss[2]]), opt_loss[1]],
    J_opt = J[opt_loss[2]],
    sigma_opt = sigma[opt_loss[1]],
    call = cl
  )

  class(out) <- "spectral"
  out
}
