% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict.lhss}
\alias{predict.lhss}
\title{Obtain predicted density ratio values from a \code{lhss} object}
\usage{
\method{predict}{lhss}(
  object,
  newdata = NULL,
  sigma = c("sigmaopt", "all"),
  lambda = c("lambdaopt", "all"),
  ...
)
}
\arguments{
\item{object}{A \code{lhss} object}

\item{newdata}{Optional \code{matrix} new data set to compute the density}

\item{sigma}{A scalar with the Gaussian kernel width}

\item{lambda}{A scalar with the regularization parameter}

\item{...}{Additional arguments to be passed to the function}
}
\value{
An array with predicted density ratio values from possibly new data,
but otherwise the numerator samples.
}
\description{
Obtain predicted density ratio values from a \code{lhss} object
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
\code{\link{predict}}, \code{\link{lhss}}
}
\keyword{lhss}
\keyword{predict}
