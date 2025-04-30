set.seed(123)
# Fit model
dr <- lhss(numerator_small, denominator_small,
           nsigma = 5, nlambda = 5, ncenters = 50)
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
