# Print a `summary.spectral` object

Print a `summary.spectral` object

## Usage

``` r
# S3 method for class 'summary.spectral'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  Object of class `summary.spectral`.

- digits:

  Number of digits to use when printing the output.

- ...:

  further arguments on how to format the number of digits.

## Value

`invisble` The inputted `summary.spectral` object.

## See also

[`print`](https://rdrr.io/r/base/print.html),
[`summary.spectral`](https://thomvolker.github.io/densityratio/reference/summary.spectral.md),
[`spectral`](https://thomvolker.github.io/densityratio/reference/spectral.md)

## Examples

``` r
set.seed(123)
# Fit model
dr <- spectral(numerator_small, denominator_small)
# Inspect model object
dr
#> 
#> Call:
#> spectral(df_numerator = numerator_small, df_denominator = denominator_small)
#> 
#> Kernel Information:
#>   Kernel type: Gaussian with L2 norm distances
#>   Number of kernels: 100
#>   sigma: num [1:10] 0.807 1.191 1.455 1.688 1.913 ...
#> 
#> 
#> Subspace dimension (J): num [1:50] 1 2 4 6 8 10 11 13 15 17 ...
#> 
#> Optimal sigma: 3.582214
#> Optimal subspace: 8
#> Optimal kernel weights (cv): num [1:8] 1.0045 -0.6689 -0.0938 0.8499 0.0228 ...
#>  
# Obtain summary of model object
summary(dr)
#> 
#> Call:
#> spectral(df_numerator = numerator_small, df_denominator = denominator_small)
#> 
#> Kernel Information:
#>   Kernel type: Gaussian with L2 norm distances
#>   Number of kernels: 100
#> 
#> Optimal sigma: 3.582214
#> Optimal subspace: 8
#> Optimal kernel weights (cv): num [1:8] 1.0045 -0.6689 -0.0938 0.8499 0.0228 ...
#>  
#> Pearson divergence between P(nu) and P(de): 0.8063
#> For a two-sample homogeneity test, use 'summary(x, test = TRUE)'.
#> 
# Plot model object
plot(dr)
#> Warning: Negative estimated density ratios for 22 observation(s) converted to 0.01 before applying logarithmic transformation
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.

# Plot density ratio for each variable individually
plot_univariate(dr)
#> Warning: Negative estimated density ratios for 22 observation(s) converted to 0.01 before applying logarithmic transformation
#> [[1]]

#> 
#> [[2]]

#> 
#> [[3]]

#> 
# Plot density ratio for each pair of variables
plot_bivariate(dr)
#> Warning: Negative estimated density ratios for 22 observation(s) converted to 0.01 before applying logarithmic transformation
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
#> [1,] 2.8779372
#> [2,] 3.8658889
#> [3,] 3.7674156
#> [4,] 4.8603842
#> [5,] 0.6970124
#> [6,] 2.1671079
#> 
# Fit model with custom parameters
spectral(numerator_small, denominator_small, sigma = 2)
#> 
#> Call:
#> spectral(df_numerator = numerator_small, df_denominator = denominator_small,     sigma = 2)
#> 
#> Kernel Information:
#>   Kernel type: Gaussian with L2 norm distances
#>   Number of kernels: 100
#>   sigma: num 2
#> 
#> 
#> Subspace dimension (J): num [1:50] 1 2 4 6 8 10 11 13 15 17 ...
#> 
#> Optimal sigma: 2
#> Optimal subspace: 4
#> Optimal kernel weights (cv): num [1:4] 0.98 -0.8324 -0.0561 0.6471
#>  
```
