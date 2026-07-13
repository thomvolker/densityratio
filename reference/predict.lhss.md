# Obtain predicted density ratio values from a `lhss` object

Obtain predicted density ratio values from a `lhss` object

## Usage

``` r
# S3 method for class 'lhss'
predict(
  object,
  newdata = NULL,
  sigma = c("sigmaopt", "all"),
  lambda = c("lambdaopt", "all"),
  ...
)
```

## Arguments

- object:

  A `lhss` object

- newdata:

  Optional `matrix` new data set to compute the density

- sigma:

  A scalar with the Gaussian kernel width

- lambda:

  A scalar with the regularization parameter

- ...:

  Additional arguments to be passed to the function

## Value

An array with predicted density ratio values from possibly new data, but
otherwise the numerator samples.

## See also

[`predict`](https://rdrr.io/r/stats/predict.html),
[`lhss`](https://thomvolker.github.io/densityratio/reference/lhss.md)

## Examples

``` r
set.seed(123)
# Fit model (minimal example to limit computation time)
dr <- lhss(numerator_small, denominator_small,
           nsigma = 5, nlambda = 3, ncenters = 50, maxit = 100)
# Inspect model object
dr
#> 
#> Call:
#> lhss(df_numerator = numerator_small, df_denominator = denominator_small,     nsigma = 5, nlambda = 3, ncenters = 50, maxit = 100)
#> 
#> Kernel Information:
#>   Kernel type: Gaussian with L2 norm distances
#>   Number of kernels: 50
#>   sigma: num [1:5, 1:3] 0.0682 0.4267 0.8044 1.3056 2.2871 ...
#> 
#> Regularization parameter (lambda): num [1:3] 1e+03 1e+00 1e-03
#> 
#> Subspace dimension (m): 1
#> Optimal sigma: 1.497984
#> Optimal lambda: 1
#> Optimal kernel weights (loocv): num [1:51] 1.8373 -0.09066 0.07367 -0.00528 0.06262 ...
#>  
# Obtain summary of model object
summary(dr)
#> 
#> Call:
#> lhss(df_numerator = numerator_small, df_denominator = denominator_small,     nsigma = 5, nlambda = 3, ncenters = 50, maxit = 100)
#> 
#> Kernel Information:
#>   Kernel type: Gaussian with L2 norm distances
#>   Number of kernels: 50
#> 
#> Subspace dimension (m): 1
#> Optimal sigma: 1.497984
#> Optimal lambda: 1
#> Optimal kernel weights (loocv): num [1:51] 1.8373 -0.09066 0.07367 -0.00528 0.06262 ...
#>  
#> Pearson divergence between P(nu) and P(de): 0.4685
#> For a two-sample homogeneity test, use 'summary(x, test = TRUE)'.
#> 
# Plot model object
plot(dr)
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.

# Plot density ratio for each variable individually
plot_univariate(dr)
#> [[1]]

#> 
#> [[2]]

#> 
#> [[3]]

#> 
# Plot density ratio for each pair of variables
plot_bivariate(dr)
#> [[1]]

#> 
#> [[2]]

#> 
#> [[3]]

#> 
# Predict density ratio and inspect first 6 predictions
head(predict(dr))
#> , , 1
#> 
#>           [,1]
#> [1,] 1.8470376
#> [2,] 2.7363298
#> [3,] 2.5131903
#> [4,] 2.7363216
#> [5,] 0.5834822
#> [6,] 1.8982260
#> 
```
