#' Least-squares heterodistributional subspace search
#'
#' @param df_numerator \code{data.frame} with exclusively numeric variables with
#' the numerator samples
#' @param df_denominator \code{data.frame} with exclusively numeric variables
#' with the denominator samples (must have the same variables as
#' \code{df_denominator})
#' @param m Scalar indicating the dimensionality of the reduced subspace
#' @param sigma \code{NULL} or a scalar value to determine the bandwidth of the
#' Gaussian kernel gram matrix. If \code{NULL}, sigma is the median Euclidean
#' interpoint distance.
#' @param lambda \code{NULL} or a scalar value to determine the regularization
#' imposed on the Gaussian kernel gram matrix of the denominator samples. If
#' \code{NULL}, \code{lambda} is chosen to be \eqn{\sqrt{N}}.
#' @param ncenters Maximum number of Gaussian centers in the kernel gram
#' matrix. Defaults to all numerator samples.
#' @param centers Numeric matrix with the same variables as \code{nu} and
#' \code{de} that are used as Gaussian centers in the kernel Gram matrix. By
#' default, the matrix \code{nu} is used as the matrix with Gaussian centers.
#' @param maxit Maximum number of iterations in the updating scheme.
#' @importFrom expm expm
#'
#' @return \code{lhss} returns \code{rhat}, the estimated density ratio.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' N <- 1000
#' X <- cbind(rnorm(N), rnorm(N, 0, 0.5))
#' Y <- cbind(rnorm(N), sample(rep(c(-1, 1), times = N/2)) + rnorm(N))
#' out <- lhss(X, Y, m = 1, ncenters = 100)
#' }


lhss <- function(df_numerator, df_denominator, m = 1, nsigma = 10,
                 sigma_quantile = NULL, sigma = NULL, nlambda = 10,
                 lambda = NULL, ncenters = 200, centers = NULL, maxit = 200,
                 parallel = FALSE, nthreads = NULL, progressbar = TRUE) {
  # TODO: Add checks

  cl <- match.call()
  nu <- as.matrix(df_numerator)
  de <- as.matrix(df_denominator)

  check.dataform(nu, de)
  centers   <- check.centers(nu, centers, ncenters)
  symmetric <- check.symmetric(nu, centers)
  parallel  <- check.parallel(parallel, nthreads, sigma, lambda)
  nthreads  <- check.threads(parallel, nthreads)
  nsigma    <- check.nsigma.lhss(nsigma, sigma, sigma_quantile)
  lambda    <- check.lambda(nlambda, lambda)

  p <- ncol(nu)
  n_nu <- nrow(nu)
  n_de <- nrow(de)

  U <- make_UV(matrix(1, p, m))[, seq_len(m), drop = FALSE]
  # old code; kept as backup
  # U_init <- matrix(1, p, m)
  # U <- .update_UV(U_init, m, p)$U

  nu_u  <- nu %*% U
  de_u  <- de %*% U
  ce_u  <- centers %*% U
  dist_nu_u <- distance(nu_u, ce_u)
  dist_de_u <- distance(de_u, ce_u)
  sigma <- check.sigma(nsigma, sigma_quantile, sigma, dist_nu_u)

  ulsif_init <- compute_ulsif(dist_nu_u, dist_de_u, sigma, lambda, FALSE, 0, FALSE)
  min_score  <- which.min(ulsif_init$loocv_score) - 1
  alpha_min  <- c(lambda = min_score %% length(lambda) + 1,
                  sigma  = min_score %/% length(lambda) + 1)
  alphah <- pmax(0, ulsif_init$alpha[, alpha_min["lambda"], alpha_min["sigma"]])

  sigma_opt <- sigma[alpha_min["sigma"]]
  alpha_opt <- alphah
  PD_opt    <- crossprod(hhat, alphah) - 0.5
  U_opt     <- U
  score_opt <- min(ulsif_init$loocv_score)

  decPDcount <- 0
  conv <- FALSE
  iter <- 0

  one_n_nu   <- rep(1, n_nu)
  one_n_de   <- rep(1, n_de)
  one_n_ce_T <- rep(1, ncenters) |> t()
  one_m_T    <- rep(1, m) |> t()

  while (!conv) {
    iter <- iter + 1
    cat(paste0("\r Iteration: ", iter))

    UV <- .update_UV(U, m, p)
    U <- UV$U
    V <- UV$V

    nu_u <- nu %*% U
    de_u <- de %*% U
    ce_u <- centers %*% U

    dPd1 <- matrix(0, p, m)

    dist_nu_u <- distance(nu_u, ce_u)
    sigma <- median_distance(dist_nu_u)
    Ktemp <- kernel_gaussian(dist_nu_u, sigma)


    for (i in 1:ncenters) {
      temp <- -(nu_u - tcrossprod(one_n_nu, ce_u[i, ])) * (Ktemp[, i] %*% one_m_T)
      temp1 <- (nu - tcrossprod(one_n_nu, centers[i, ]))
      h <- crossprod(temp, temp1)

      dPd1 <- dPd1 + t(h) * alphah[i]
    }

    dPd1 <- dPd1 / n_nu / sigma

    Ktemp <- kernel_gaussian(distance(de_u, ce_u), sigma)
    dPd2 <- matrix(0, p, m)

    for (i in 1:ncenters) {
      temp11 <- (de_u - tcrossprod(one_n_de, ce_u[i, ])) * (Ktemp[, i] %*% one_m_T)
      tempx1 <- (de - tcrossprod(one_n_de, centers[i, ]))

      for (j in 1:m) {
        T11 <- temp11[, j] %*% one_n_ce_T * Ktemp
        dPd2[, j] <- dPd2[, j] - crossprod(tempx1, T11) %*% alphah * alphah[i] * 2
      }
    }

    dPd2 <- dPd2 / n_de / sigma

    dPd <- t(dPd1 - dPd2/2)

    dM <- rbind(
      cbind(matrix(0, m, m), dPd %*% V),
      cbind(-t(dPd %*% V), matrix(0, p - m, p - m))
    )

    M  <- dM * 1/maxit * (maxit - (iter - 1))
    eM <- expm::expm(M)
    U  <- cbind(diag(m), matrix(0, m, p - m)) %*% eM %*% rbind(t(U), t(V))
    U  <- t(U)

    nu_u <- nu %*% U
    de_u <- de %*% U
    ce_u <- centers %*% U

    dist_nu_u <- distance(nu_u, ce_u)
    sigma <- median_distance(dist_nu_u)

    phi_nu <- kernel_gaussian(dist_nu_u, sigma)
    phi_de <- kernel_gaussian(distance(de_u, ce_u), sigma)
    Hhat   <- crossprod(phi_de) / n_de
    hhat   <- colMeans(phi_nu)

    alphat <- compute_ulsif(Hhat, hhat, lambda, FALSE, 0)
    alphah <- pmax(0, alphat)
    PD <- crossprod(hhat, alphah) - 0.5

    # TODO: investigate better stopping criterions (e.g., dynamic stepsize update (i.e., half current M))
    if (PD > PD_opt) {
      U_opt <- U
      alphah_opt <- alphah
      PD_opt <- PD
      sigmatopt <- sigma
      decPDcount <- max(0, decPDcount - 1)
    } else {
      decPDcount <- decPDcount + 1
    }

    if (iter == maxit | decPDcount > 20) {
      conv <- TRUE
    }
  }

  list(U = U_opt,
       alpha = alphah_opt,
       sigma = sigmatopt,
       PD = PD_opt,
       call = cl)
}

.update_UV <- function(U, m, p) {
  QR <- qr(U)
  Q <- qr.Q(QR, complete = TRUE) * sign(qr.R(QR, complete = TRUE)[1])
  U <- Q[, seq_len(m), drop = FALSE]
  V <- Q[, seq(m+1, p), drop = FALSE]
  list(U = U, V = V)
}
