
\dontshow{
if (requireNamespace("RcppArmadillo", quietly = TRUE)) {
  RcppArmadillo::armadillo_throttle_cores(1)
}
}

set.seed(123)
# Fit model
dr <- rulsif(numerator_small, denominator_small)
# Inspect model object
dr
# With identical Gaussian centers, the following models are identical
rulsif(numerator_small, denominator_small, a = 0)
ulsif(numerator_small, denominator_small)

# Additional helper functions
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
rulsif(numerator_small, denominator_small, sigma = 2, lambda = 2)


\dontshow{
if(requireNamespace("RcppArmadillo", quietly = TRUE)) {
RcppArmadillo::armadillo_reset_cores()
}
}