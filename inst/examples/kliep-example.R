set.seed(123)
# Fit model
dr <- kliep(numerator_small, denominator_small)
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
kliep(numerator_small, denominator_small,
      nsigma = 20, ncenters = 100, nfold = 20,
      epsilon = 10^{3:-5}, maxit = 1000)
