#' Least-squares heterodistributional subspace search (naive implementation)
#'
#' @param nu Numeric matrix with numerator samples
#' @param de Numeric matrix with denominator samples (must have the same
#' variables as \code{nu})
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
#' @importFrom pracma householder
#' @importFrom expm expm
#'
#' @return \code{lhss_naive} returns \code{rhat}, the estimated density ratio.
#'
#' @examples
#' set.seed(1)
#' N <- 1000
#' X <- cbind(rnorm(N), rnorm(N, 0, 0.5))
#' Y <- cbind(rnorm(N), sample(rep(c(-1, 1), times = N/2)) + rnorm(N))
#' out <- lhss_naive(X, Y, m = 1, ncenters = 100)


lhss_naive <- function(nu, de, m = 1, sigma = NULL, lambda = 1,
                       ncenters = nrow(nu), centers = NULL, maxit = 200) {
  # TO DO: Add checks

  nu <- as.matrix(nu)
  de <- as.matrix(de)

  p <- ncol(nu)
  n_nu <- nrow(nu)
  n_de <- nrow(de)

  U <- matrix(1, p, m)
  QR <- pracma::householder(U)
  Q <- QR$Q * sign(QR$R[1])
  U <- Q[ , 1:m, drop = F]

  if (is.null(centers)) {
    if (ncenters < nrow(nu)) {
      centers <- nu[sample(n_nu, ncenters), ]
    } else {
      centers <- nu
    }
  } else {
    centers <- as.matrix(centers)
    ncenters <- nrow(centers)
    if (!is.numeric(centers) | ! p == ncol(centers)) {
      stop("If centers are provided, they must have the same variables as the numerator samples")
    }
  }

  nu_u  <- nu %*% U
  de_u  <- de %*% U
  ce_u  <- centers %*% U
  dist_nu_u <- distance(nu_u, ce_u)
  sigma <- median_distance(dist_nu_u)

  phi_nu <- kernel_gaussian(dist_nu_u, sigma)
  phi_de <- kernel_gaussian(distance(de_u, ce_u), sigma)
  Hhat   <- crossprod(phi_de) / n_de
  hhat   <- colMeans(phi_nu)

  eyeH   <- lambda * diag(ncenters)
  alphat <- solve(Hhat + eyeH) %*% hhat
  alphah <- pmax(0, alphat)
  PD_opt <- crossprod(hhat, alphah) - 0.5
  U_opt  <- U

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

    QR <- pracma::householder(U)
    Q  <- QR$Q * sign(QR$R[1])
    U  <- Q[, 1:m, drop = FALSE]
    V  <- Q[, (m + 1):p, drop = FALSE]

    nu_u <- nu %*% U
    de_u <- de %*% U
    ce_u <- centers %*% U
    h <- rep(0, p)

    dPd1 <- matrix(0, p, m)


    dist_nu_u <- distance(nu_u, ce_u)
    sigma <- median_distance(dist_nu_u)
    Ktemp <- kernel_gaussian(dist_nu_u, sigma)


    for (i in 1:ncenters) {
      temp <- -(nu_u - one_n_nu %*% t(ce_u[i, ])) * (Ktemp[, i] %*% one_m_T)
      temp1 <- (nu - one_n_nu %*% t(centers[i, ]))
      h <- crossprod(temp, temp1)

      dPd1 <- dPd1 + t(h) * alphah[i]
    }

    dPd1 <- dPd1 / n_nu / sigma

    Ktemp <- kernel_gaussian(distance(de_u, ce_u), sigma)
    dPd2 <- matrix(0, p, m)

    for (i in 1:ncenters) {
      temp11 <- (de_u - one_n_de %*% t(ce_u[i, ])) * (Ktemp[, i] %*% one_m_T)
      tempx1 <- (de - one_n_de %*% t(centers[i, ]))

      for (j in 1:m) {
        T11 <- temp11[, j] %*% one_n_ce_T * Ktemp
        dPd2[, j] <- dPd2[, j] - t(tempx1) %*% T11 %*% alphah * alphah[i] * 2
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

    alphat <- solve(Hhat + eyeH) %*% hhat
    alphah <- pmax(0, alphat)
    PD <- crossprod(hhat, alphah) - 0.5

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
       PD = PD_opt)
}
