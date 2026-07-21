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
           nsigma = 3, lambda = c(0.1, 1), ncenters = 50, maxit = 100)
# Inspect model object
dr
#> 
#> Call:
#> lhss(df_numerator = numerator_small, df_denominator = denominator_small,     nsigma = 3, lambda = c(0.1, 1), ncenters = 50, maxit = 100)
#> 
#> Kernel Information:
#>   Kernel type: Gaussian with L2 norm distances
#>   Number of kernels: 50
#>   sigma: num [1:3, 1:2] 0.084 0.9371 2.6764 0.0703 0.9633 ...
#> 
#> Regularization parameter (lambda): num [1:2] 0.1 1
#> 
#> Subspace dimension (m): 1
#> Optimal sigma: 2.676351
#> Optimal lambda: 0.1
#> Optimal kernel weights (loocv): num [1:51] 6.6612 -0.308 -0.0378 -0.0356 0.0108 ...
#>  
# Obtain summary of model object
summary(dr)
#> 
#> Call:
#> lhss(df_numerator = numerator_small, df_denominator = denominator_small,     nsigma = 3, lambda = c(0.1, 1), ncenters = 50, maxit = 100)
#> 
#> Kernel Information:
#>   Kernel type: Gaussian with L2 norm distances
#>   Number of kernels: 50
#> 
#> Subspace dimension (m): 1
#> Optimal sigma: 2.676351
#> Optimal lambda: 0.1
#> Optimal kernel weights (loocv): num [1:51] 6.6612 -0.308 -0.0378 -0.0356 0.0108 ...
#>  
#> Pearson divergence between P(nu) and P(de): 0.6365
#> For a two-sample homogeneity test, use 'summary(x, test = TRUE)'.
#> 
# Plot model object
plot(dr)
#> Warning: Negative estimated density ratios for 13 observation(s) converted to 0.01 before applying logarithmic transformation
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.

# Plot density ratio for each variable individually
plot_univariate(dr)
#> Warning: Negative estimated density ratios for 13 observation(s) converted to 0.01 before applying logarithmic transformation
#> [[1]]

#> 
#> [[2]]

#> 
#> [[3]]

#> 
# Plot density ratio for each pair of variables
plot_bivariate(dr)
#> Warning: Negative estimated density ratios for 13 observation(s) converted to 0.01 before applying logarithmic transformation
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
#> [1,] 2.1662109
#> [2,] 3.6734863
#> [3,] 2.8023392
#> [4,] 4.1434823
#> [5,] 0.1595074
#> [6,] 1.2642885
#> 
```
