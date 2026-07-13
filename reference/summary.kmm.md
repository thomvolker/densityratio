# Extract summary from `kmm` object, including two-sample significance test for homogeneity of the numerator and denominator samples

Extract summary from `kmm` object, including two-sample significance
test for homogeneity of the numerator and denominator samples

## Usage

``` r
# S3 method for class 'kmm'
summary(
  object,
  test = FALSE,
  n_perm = 100,
  parallel = FALSE,
  cluster = NULL,
  min_pred = 1e-06,
  ...
)
```

## Arguments

- object:

  Object of class `kmm`

- test:

  logical indicating whether to statistically test for homogeneity of
  the numerator and denominator samples.

- n_perm:

  Scalar indicating number of permutation samples

- parallel:

  `logical` indicating to run the permutation test in parallel

- cluster:

  `NULL` or a cluster object created by `makeCluster`. If `NULL` and
  `parallel = TRUE`, it uses the number of available cores minus 1.

- min_pred:

  Scalar indicating the minimum value for the predicted density ratio
  values (used in the divergence statistic) to avoid negative density
  ratio values.

- ...:

  further arguments passed to or from other methods.

## Value

Summary of the fitted density ratio model

## Examples

``` r
set.seed(123)
# Fit model
dr <- kmm(numerator_small, denominator_small)
# Inspect model object
dr
#> 
#> Call:
#> kmm(df_numerator = numerator_small, df_denominator = denominator_small)
#> 
#> Kernel Information:
#>   Kernel type: Gaussian with L2 norm distances
#>   Number of kernels: 150
#>   sigma: num [1:10] 0.801 1.2 1.483 1.723 1.954 ...
#> 
#> Optimal sigma (5-fold cv): 3.67
#> Optimal kernel weights (5-fold cv):  num [1:150, 1] 0.23 0.416 -0.166 1.512 0.831 ...
#> 
#> Optimization parameters:
#>   Optimization method:  Unconstrained 
#> 
# Obtain summary of model object
summary(dr)
#> 
#> Call:
#> kmm(df_numerator = numerator_small, df_denominator = denominator_small)
#> 
#> Kernel Information:
#>   Kernel type: Gaussian with L2 norm distances
#>   Number of kernels: 150
#> Optimal sigma: 3.669758
#> Optimal kernel weights: num [1:150, 1] 0.23 0.416 -0.166 1.512 0.831 ...
#>  
#> Pearson divergence between P(nu) and P(de): 0.9439
#> For a two-sample homogeneity test, use 'summary(x, test = TRUE)'.
#> 
# Plot model object
plot(dr)
#> Warning: Negative estimated density ratios for 19 observation(s) converted to 0.01 before applying logarithmic transformation
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.

# Plot density ratio for each variable individually
plot_univariate(dr)
#> Warning: Negative estimated density ratios for 19 observation(s) converted to 0.01 before applying logarithmic transformation
#> [[1]]

#> 
#> [[2]]

#> 
#> [[3]]

#> 
# Plot density ratio for each pair of variables
plot_bivariate(dr)
#> Warning: Negative estimated density ratios for 19 observation(s) converted to 0.01 before applying logarithmic transformation
#> [[1]]

#> 
#> [[2]]

#> 
#> [[3]]

#> 
# Predict density ratio and inspect first 6 predictions
head(predict(dr))
#>           [,1]
#> [1,] 3.1261579
#> [2,] 4.0233887
#> [3,] 3.6868339
#> [4,] 5.5934888
#> [5,] 0.6302996
#> [6,] 1.5225886
# Fit model with custom parameters
kmm(numerator_small, denominator_small,
    nsigma = 5, ncenters = 100, nfold = 10,
    constrained = TRUE)
#> 
#> Call:
#> kmm(df_numerator = numerator_small, df_denominator = denominator_small,     constrained = TRUE, nsigma = 5, ncenters = 100, nfold = 10)
#> 
#> Kernel Information:
#>   Kernel type: Gaussian with L2 norm distances
#>   Number of kernels: 100
#>   sigma: num [1:5] 0.811 1.577 2.094 2.66 3.706
#> 
#> Optimal sigma (10-fold cv): 2.094
#> Optimal kernel weights (10-fold cv):  num [1:100, 1] -0.000498 -0.000999 -0.001187 -0.001022 -0.000275 ...
#> 
#> Optimization parameters:
#>   Optimization method:  Constrained 
#> 
```
