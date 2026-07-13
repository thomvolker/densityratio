# Obtain predicted density ratio values from a `kmm` object

Obtain predicted density ratio values from a `kmm` object

## Usage

``` r
# S3 method for class 'kmm'
predict(object, newdata = NULL, sigma = c("sigmaopt", "all"), ...)
```

## Arguments

- object:

  A `kmm` object

- newdata:

  Optional `matrix` new data set to compute the density

- sigma:

  A scalar with the Gaussian kernel width

- ...:

  Additional arguments to be passed to the function

## Value

An array with predicted density ratio values from possibly new data, but
otherwise the numerator samples.

## See also

[`predict`](https://rdrr.io/r/stats/predict.html),
[`kmm`](https://thomvolker.github.io/densityratio/reference/kmm.md)

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
