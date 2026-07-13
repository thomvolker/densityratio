# Print a `summary.lhss` object

Print a `summary.lhss` object

## Usage

``` r
# S3 method for class 'summary.lhss'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  Object of class `summary.lhss`.

- digits:

  Number of digits to use when printing the output.

- ...:

  further arguments on how to format the number of digits.

## Value

`invisble` The inputted `summary.lhss` object.

## See also

[`print`](https://rdrr.io/r/base/print.html),
[`summary.lhss`](https://thomvolker.github.io/densityratio/reference/summary.lhss.md),
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
