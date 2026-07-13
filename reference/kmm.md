# Kernel mean matching approach to density ratio estimation

Kernel mean matching approach to density ratio estimation

## Usage

``` r
kmm(
  df_numerator,
  df_denominator,
  scale = "numerator",
  constrained = FALSE,
  nsigma = 10,
  sigma_quantile = NULL,
  sigma = NULL,
  ncenters = 200,
  centers = NULL,
  cv = TRUE,
  nfold = 5,
  parallel = FALSE,
  nthreads = NULL,
  progressbar = TRUE,
  osqp_settings = NULL,
  cluster = NULL
)
```

## Arguments

- df_numerator:

  `data.frame` with exclusively numeric variables with the numerator
  samples

- df_denominator:

  `data.frame` with exclusively numeric variables with the denominator
  samples (must have the same variables as `df_denominator`)

- scale:

  `"numerator"`, `"denominator"`, or `NULL`, indicating whether to
  standardize each numeric variable according to the numerator means and
  standard deviations, the denominator means and standard deviations, or
  apply no standardization at all.

- constrained:

  `logical` equals `FALSE` to use unconstrained optimization, `TRUE` to
  use constrained optimization. Defaults to `FALSE`.

- nsigma:

  Integer indicating the number of sigma values (bandwidth parameter of
  the Gaussian kernel gram matrix) to use in cross-validation.

- sigma_quantile:

  `NULL` or numeric vector with probabilities to calculate the quantiles
  of the distance matrix to obtain sigma values. If `NULL`, `nsigma`
  values between `0.25` and `0.75` are used.

- sigma:

  `NULL` or a scalar value to determine the bandwidth of the Gaussian
  kernel gram matrix. If `NULL`, `nsigma` values between `0.25` and
  `0.75` are used.

- ncenters:

  `integer` Maximum number of Gaussian centers in the kernel gram
  matrix.

- centers:

  `NULL` or `data.frame` with the same dimensions as the data,
  indicating the centers for the Gaussian kernel gram matrix.

- cv:

  Logical indicating whether or not to do cross-validation

- nfold:

  Number of cross-validation folds used in order to calculate the
  optimal `sigma` value (default is 5-fold cv).

- parallel:

  logical indicating whether to use parallel processing in the
  cross-validation scheme.

- nthreads:

  `NULL` or integer indicating the number of threads to use for parallel
  processing. If parallel processing is enabled, it defaults to the
  number of available threads minus one.

- progressbar:

  Logical indicating whether or not to display a progressbar.

- osqp_settings:

  Optional: settings to pass to the `osqp` solver for constrained
  optimization.

- cluster:

  Optional: a cluster object to use for parallel processing, see
  [`parallel::makeCluster`](https://rdrr.io/r/parallel/makeCluster.html).

## Value

`kmm`-object, containing all information to calculate the density ratio
using optimal sigma and optimal weights.

## References

Huang, J., Smola, A. J., Gretton, A., Borgwardt, K. M., & Schölkopf, B.
(2006). Correcting sample selection bias by unlabeled data. In *Advances
in Neural Information Processing Systems*, edited by B. Schölkopf, J.
Platt and T. Hoffman. Available from
<https://proceedings.neurips.cc/paper/2006/hash/a2186aa7c086b46ad4e8bf81e2a3a19b-Abstract.html>.

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
