# A histogram of density ratio estimates

Creates a histogram of the density ratio estimates. Useful to understand
the distribution of estimated density ratios in each sample, or compare
it among samples. It is the default plotting method for density ratio
objects.

## Usage

``` r
dr.histogram(
  x,
  samples = "both",
  logscale = TRUE,
  binwidth = NULL,
  bins = NULL,
  tol = 0.01,
  ...
)

# S3 method for class 'ulsif'
plot(
  x,
  samples = "both",
  logscale = TRUE,
  binwidth = NULL,
  bins = NULL,
  tol = 0.01,
  ...
)

# S3 method for class 'kliep'
plot(
  x,
  samples = "both",
  logscale = TRUE,
  binwidth = NULL,
  bins = NULL,
  tol = 0.01,
  ...
)

# S3 method for class 'kmm'
plot(
  x,
  samples = "both",
  logscale = TRUE,
  binwidth = NULL,
  bins = NULL,
  tol = 0.01,
  ...
)

# S3 method for class 'spectral'
plot(
  x,
  samples = "both",
  logscale = TRUE,
  binwidth = NULL,
  bins = NULL,
  tol = 0.01,
  ...
)

# S3 method for class 'lhss'
plot(
  x,
  samples = "both",
  logscale = TRUE,
  binwidth = NULL,
  bins = NULL,
  tol = 0.01,
  ...
)

# S3 method for class 'naivedensityratio'
plot(
  x,
  samples = "both",
  logscale = TRUE,
  binwidth = NULL,
  bins = NULL,
  tol = 0.01,
  ...
)
```

## Arguments

- x:

  Density ratio object created with e.g.,
  [`kliep()`](https://thomvolker.github.io/densityratio/reference/kliep.md),
  [`ulsif()`](https://thomvolker.github.io/densityratio/reference/ulsif.md),
  or
  [`naive()`](https://thomvolker.github.io/densityratio/reference/naive.md)

- samples:

  Character string indicating whether to plot the 'numerator',
  'denominator', or 'both' samples. Default is 'both'.

- logscale:

  Logical indicating whether to plot the density ratio estimates on a
  log scale. Default is TRUE.

- binwidth:

  Numeric indicating the width of the bins, passed on to `ggplot2`.

- bins:

  Numeric indicating the number of bins. Overriden by binwidth, and
  passed on to `ggplot2`.

- tol:

  Numeric indicating the tolerance: values below this value will be set
  to the tolerance value, for legibility of the plots

- ...:

  Additional arguments passed on to
  [`predict()`](https://rdrr.io/r/stats/predict.html).

## Value

A histogram of density ratio estimates.

A histogram of density ratio estimates.

A histogram of density ratio estimates.

A histogram of density ratio estimates.

A histogram of density ratio estimates.

A histogram of density ratio estimates.

A histogram of density ratio estimates.

## See also

[`ulsif`](https://thomvolker.github.io/densityratio/reference/ulsif.md)
for example usage

[`kliep`](https://thomvolker.github.io/densityratio/reference/kliep.md)
for example usage

[`kmm`](https://thomvolker.github.io/densityratio/reference/kmm.md) for
example usage

[`spectral`](https://thomvolker.github.io/densityratio/reference/spectral.md)
for example usage

[`lhss`](https://thomvolker.github.io/densityratio/reference/lhss.md)
for example usage

[`naive`](https://thomvolker.github.io/densityratio/reference/naive.md)
for example usage
