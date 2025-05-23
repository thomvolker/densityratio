% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lhss.R
\name{lhss}
\alias{lhss}
\title{Least-squares heterodistributional subspace search}
\usage{
lhss(
  df_numerator,
  df_denominator,
  m = NULL,
  intercept = TRUE,
  scale = "numerator",
  nsigma = 10,
  sigma_quantile = NULL,
  sigma = NULL,
  nlambda = 10,
  lambda = NULL,
  ncenters = 200,
  centers = NULL,
  maxit = 200,
  progressbar = TRUE
)
}
\arguments{
\item{df_numerator}{\code{data.frame} with exclusively numeric variables with
the numerator samples}

\item{df_denominator}{\code{data.frame} with exclusively numeric variables
with the denominator samples (must have the same variables as
\code{df_denominator})}

\item{m}{Scalar indicating the dimensionality of the reduced subspace}

\item{intercept}{\code{logical} Indicating whether to include an intercept
term in the model. Defaults to \code{TRUE}.}

\item{scale}{\code{"numerator"}, \code{"denominator"}, or \code{NULL},
indicating whether to standardize each numeric variable according to the
numerator means and standard deviations, the denominator means and standard
deviations, or apply no standardization at all.}

\item{nsigma}{Integer indicating the number of sigma values (bandwidth
parameter of the Gaussian kernel gram matrix) to use in cross-validation.}

\item{sigma_quantile}{\code{NULL} or numeric vector with probabilities to
calculate the quantiles of the distance matrix to obtain sigma values. If
\code{NULL}, \code{nsigma} values between \code{0.05} and \code{0.95} are
used.}

\item{sigma}{\code{NULL} or a scalar value to determine the bandwidth of the
Gaussian kernel gram matrix. If \code{NULL}, \code{nsigma} values between
\code{0.05} and \code{0.95} are used.}

\item{nlambda}{Integer indicating the number of \code{lambda} values
(regularization parameter), by default, \code{lambda} is set to
\code{10^seq(3, -3, length.out = nlambda)}.}

\item{lambda}{\code{NULL} or numeric vector indicating the lambda values to
use in cross-validation}

\item{ncenters}{Maximum number of Gaussian centers in the kernel gram
matrix. Defaults to all numerator samples.}

\item{centers}{Numeric matrix with the same variables as \code{nu} and
\code{de} that are used as Gaussian centers in the kernel Gram matrix. By
default, the matrix \code{nu} is used as the matrix with Gaussian centers.}

\item{maxit}{Maximum number of iterations in the updating scheme.}

\item{progressbar}{Logical indicating whether or not to display a progressbar.}
}
\value{
\code{lhss}-object, containing all information to calculate the
density ratio using optimal sigma, optimal lambda and optimal weights.
}
\description{
Least-squares heterodistributional subspace search
}
\examples{
set.seed(123)
# Fit model
dr <- naive(numerator_small, denominator_small)
# Inspect model object
dr
# Obtain summary of model object
summary(dr)
# Plot model object
plot(dr)
# Plot density ratio for each variable individually
plot_univariate(dr)
# Plot density ratio for each pair of variables
plot_bivariate(dr)
# Predict density ratio and inspect first 6 predictions
head(predict(dr))
# Fit model with custom parameters
naive(numerator_small, denominator_small, m=2, kernel="epanechnikov")
}
\references{
Sugiyama, M., Yamada, M., Von Bünau, P., Suzuki, T., Kanamori, T.
& Kawanabe, M. (2011). Direct density-ratio estimation with dimensionality
reduction via least-squares hetero-distributional subspace search. \if{html}{\out{<i>}} Neural
Networks\if{html}{\out{</i>}}, 24, 183-198. \doi{10.1016/j.neunet.2010.10.005}.
}
