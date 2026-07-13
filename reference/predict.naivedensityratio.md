# Obtain predicted density ratio values from a `naivedensityratio` object

Obtain predicted density ratio values from a `naivedensityratio` object

## Usage

``` r
# S3 method for class 'naivedensityratio'
predict(object, newdata = NULL, log = FALSE, tol = 1e-06, ...)
```

## Arguments

- object:

  A `naive` object

- newdata:

  Optional `matrix` new data set to compute the density

- log:

  A logical indicating whether to return the log of the density ratio

- tol:

  Minimal density value to avoid numerical issues

- ...:

  Additional arguments to be passed to the function

## Value

An array with predicted density ratio values from possibly new data, but
otherwise the numerator samples.

## See also

[`predict`](https://rdrr.io/r/stats/predict.html),
[`naive`](https://thomvolker.github.io/densityratio/reference/naive.md)

## Examples

``` r
set.seed(123)
# Fit model
dr <- naive(numerator_small, denominator_small)
# Inspect model object
dr
#> 
#> Call:
#> naive(df_numerator = numerator_small, df_denominator = denominator_small)
#> 
#> Naive density ratio
#>   Number of variables: 3
#>   Number of numerator samples: 50
#>   Number of denominator samples: 100
#>   Numerator density: num [1:50] 1.41 5.74 1.87 4.13 1.67 ...
#>   Denominator density: num [1:100] 2.93 0.071 1.065 1.59 2.115 ...
#>  
# Obtain summary of model object
summary(dr)
#> 
#> Call:
#> naive(df_numerator = numerator_small, df_denominator = denominator_small)
#> 
#> Naive density ratio estimate:
#>   Number of variables: 
#>   Number of numerator samples: 50
#>   Number of denominator samples: 100
#>   Density ratio for numerator samples: num [1:50] 0.344 1.747 0.628 1.419 0.511 ...
#>   Density ratio for denominator samples: num [1:100] 1.0751 -2.6454 0.0626 0.464 0.7493 ...
#>  
#> 
#> Squared average log density ratio difference for numerator and denominator samples (SALDRD): 13.56
#> For a two-sample homogeneity test, use 'summary(x, test = TRUE)'.
#> 
# Plot model object
plot(dr)
#> Warning: Negative estimated density ratios for 25 observation(s) converted to 0.01 before applying logarithmic transformation
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.

# Plot density ratio for each variable individually
plot_univariate(dr)
#> Warning: Negative estimated density ratios for 25 observation(s) converted to 0.01 before applying logarithmic transformation
#> [[1]]

#> 
#> [[2]]

#> 
#> [[3]]

#> 
# Plot density ratio for each pair of variables
plot_bivariate(dr)
#> Warning: Negative estimated density ratios for 25 observation(s) converted to 0.01 before applying logarithmic transformation
#> [[1]]

#> 
#> [[2]]

#> 
#> [[3]]

#> 
# Predict density ratio and inspect first 6 predictions
head(predict(dr))
#> [1] 1.410607 5.739287 1.874031 4.131255 1.666760 4.095855
# Fit model with custom parameters
naive(numerator_small, denominator_small, m=2, kernel="epanechnikov")
#> 
#> Call:
#> naive(df_numerator = numerator_small, df_denominator = denominator_small,     m = 2, kernel = "epanechnikov")
#> 
#> Naive density ratio
#>   Number of variables: 3
#>   Number of numerator samples: 50
#>   Number of denominator samples: 100
#>   Numerator density: num [1:50] 0.572 1.421 0.945 1.058 0.936 ...
#>   Denominator density: num [1:100] 1.391 1.459 0.572 0.943 1.314 ...
#>  
```
