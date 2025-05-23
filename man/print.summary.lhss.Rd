% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{print.summary.lhss}
\alias{print.summary.lhss}
\title{Print a \code{summary.lhss} object}
\usage{
\method{print}{summary.lhss}(x, digits = max(3L, getOption("digits") - 3L), ...)
}
\arguments{
\item{x}{Object of class \code{summary.lhss}.}

\item{digits}{Number of digits to use when printing the output.}

\item{...}{further arguments on how to format the number of digits.}
}
\value{
\code{invisble} The inputted \code{summary.lhss} object.
}
\description{
Print a \code{summary.lhss} object
}
\examples{
set.seed(123)
# Fit model (minimal example to limit computation time)
dr <- lhss(numerator_small, denominator_small,
           nsigma = 5, nlambda = 3, ncenters = 50, maxit = 100)
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
}
\seealso{
\code{\link{print}}, \code{\link{summary.lhss}}, \code{\link{lhss}}
}
