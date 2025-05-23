% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/naive.R
\name{naive}
\alias{naive}
\title{Naive density ratio estimation}
\usage{
naive(
  df_numerator,
  df_denominator,
  m = NULL,
  bw = "SJ",
  kernel = "gaussian",
  n = 2L^11,
  ...
)
}
\arguments{
\item{df_numerator}{\code{data.frame} with exclusively numeric variables with
the numerator samples}

\item{df_denominator}{\code{data.frame} with exclusively numeric variables
with the denominator samples (must have the same variables as
\code{df_denominator})}

\item{m}{\code{integer} Optional parameter to reduce the dimensionality of
the data in multivariate density ratio estimation problems. If missing,
the number of variables in the data is used. If set to an integer value
smaller than the number of variables, the first \code{m} principal
components are used to estimate the density ratio. If set to \code{NULL},
the square root of the number of variables is used (for consistency with
other methods).}

\item{bw}{the smoothing bandwidth to be used. See \link[stats:density]{stats::density} for more
information.}

\item{kernel}{the kernel to be used. See \link[stats:density]{stats::density} for more
information.}

\item{n}{\code{integer} the number of equally spaced points at which the
density is to be estimated. When n > 512, it is rounded up to a power of 2
during the calculations (as fast Fourier transform is used) and the final
result is interpolated by \link[stats:approxfun]{stats::approx}. So it makes sense to specify n as
a power' of two.}

\item{...}{further arguments passed to \link[stats:density]{stats::density}}
}
\value{
\code{naivedensityratio} object
}
\description{
The naive approach creates separate kernel density estimates for
the numerator and the denominator samples, and then evaluates their
ratio for the denominator samples. For multivariate data, the density ratio
is computed after a orthogonal linear transformation, such that the new
variables can be treated as independent. To reduce the dimensionality of
the PCA solution, one can set the number of components by setting the
\code{m} parameter to an integer value smaller than the number of variables.
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
\seealso{
\code{\link[stats:density]{stats::density()}}
}
