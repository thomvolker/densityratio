# Obtain predicted density ratio values from a `spectral` object

Obtain predicted density ratio values from a `spectral` object

## Usage

``` r
# S3 method for class 'spectral'
predict(
  object,
  newdata = NULL,
  sigma = c("sigmaopt", "all"),
  m = c("opt", "all"),
  ...
)
```

## Arguments

- object:

  A `spectral` object

- newdata:

  Optional `matrix` new data set to compute the density

- sigma:

  A scalar with the Gaussian kernel width

- m:

  integer indicating the dimension of the eigenvector expansion

- ...:

  Additional arguments to be passed to the function

## Value

An array with predicted density ratio values from possibly new data, but
otherwise the numerator samples.

## See also

[`predict`](https://rdrr.io/r/stats/predict.html),
[`spectral`](https://thomvolker.github.io/densityratio/reference/spectral.md)
