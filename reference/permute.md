# Single permutation

Single permutation

Single permutation statistic of `ulsif` object

Single permutation statistic of `kliep` object

Single permutation statistic of `kmm` object

Single permutation statistic of `lhss` object

Single permutation statistic of `spectral` object

Single permutation statistic of `naivedensityratio` object

## Usage

``` r
permute(object, ...)

# S3 method for class 'ulsif'
permute(object, stacked, nnu, nde, ...)

# S3 method for class 'kliep'
permute(object, stacked, nnu, nde, min_pred = sqrt(.Machine$double.eps), ...)

# S3 method for class 'kmm'
permute(object, stacked, nnu, nde, ...)

# S3 method for class 'lhss'
permute(object, stacked, nnu, nde, ...)

# S3 method for class 'spectral'
permute(object, stacked, nnu, nde, ...)

# S3 method for class 'naivedensityratio'
permute(object, stacked, nnu, nde, min_pred, max_pred, ...)
```

## Arguments

- object:

  `naivedensityratio` object

- ...:

  Additional arguments to pass through to specific permute functions.

- stacked:

  `matrix` with stacked numerator and denominator samples

- nnu:

  Scalar with numerator sample size

- nde:

  Scalar with denominator sample size

- min_pred:

  Minimum value of the predicted density ratio

- max_pred:

  Maximum value of the predicted density ratio

## Value

permutation statistic for a single permutation of the data

permutation statistic for a single permutation of the data

permutation statistic for a single permutation of the data

permutation statistic for a single permutation of the data

permutation statistic for a single permutation of the data

permutation statistic for a single permutation of the data

permutation statistic for a single permutation of the data
