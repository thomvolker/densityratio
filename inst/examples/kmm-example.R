set.seed(123)
# Fit model
dr <- kmm(numerator_small, denominator_small)
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
kmm(numerator_small, denominator_small,
    nsigma = 20, ncenters = 100, nfold = 20,
    constrained = TRUE, parallel = TRUE, nthreads = 2)
