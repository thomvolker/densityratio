
\dontshow{
if (requireNamespace("RcppArmadillo", quietly = TRUE)) {
  RcppArmadillo::armadillo_throttle_cores(1)
}
}

set.seed(123)
# Fit model (minimal example to limit computation time)
dr <- lhss(numerator_small, denominator_small,
           nsigma = 3, lambda = c(0.1, 1), ncenters = 50, maxit = 100)
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

\dontshow{
if(requireNamespace("RcppArmadillo", quietly = TRUE)) {
RcppArmadillo::armadillo_reset_cores()
}
}