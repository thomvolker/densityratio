% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{print.spectral}
\alias{print.spectral}
\title{Print a \code{spectral} object}
\usage{
\method{print}{spectral}(x, digits = max(3L, getOption("digits") - 3L), ...)
}
\arguments{
\item{x}{Object of class \code{spectral}.}

\item{digits}{Number of digits to use when printing the output.}

\item{...}{further arguments on how to format the number of digits.}
}
\value{
\code{invisble} The inputted \code{spectral} object.
}
\description{
Print a \code{spectral} object
}
\examples{
set.seed(123)
# Fit model
dr <- spectral(numerator_small, denominator_small)
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
spectral(numerator_small, denominator_small, sigma = 2)
}
\seealso{
\code{\link{print}}, \code{\link{spectral}}
}
