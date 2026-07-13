# Print a `lhss` object

Print a `lhss` object

## Usage

``` r
# S3 method for class 'lhss'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  Object of class `lhss`.

- digits:

  Number of digits to use when printing the output.

- ...:

  further arguments on how to format the number of digits.

## Value

`invisble` The inputted `lhss` object.

## See also

[`print`](https://rdrr.io/r/base/print.html),
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
#>   sigma: num [1:5, 1:3] 0.0682 0.4267 0.8044 1.3056 2.2869 ...
#> 
#> Regularization parameter (lambda): num [1:3] 1e+03 1e+00 1e-03
#> 
#> Subspace dimension (m): 1
#> Optimal sigma: 1.494882
#> Optimal lambda: 1
#> Optimal kernel weights (loocv): num [1:51] 1.82362 -0.09005 0.07498 -0.00553 0.06341 ...
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
#> Optimal sigma: 1.494882
#> Optimal lambda: 1
#> Optimal kernel weights (loocv): num [1:51] 1.82362 -0.09005 0.07498 -0.00553 0.06341 ...
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
#> [1,] 1.8378605
#> [2,] 2.7285745
#> [3,] 2.5115771
#> [4,] 2.7286374
#> [5,] 0.5922951
#> [6,] 1.9138572
#> 
```
