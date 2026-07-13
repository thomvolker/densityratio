# Changelog

## densityratio (development version)

- update loocv procedure in ulsif to work with unconstrained kernel
  weights
- update lhss to perform unconstrained optimization (aligned with loocv
  procedure of ulsif)
- small formal check to stop when negative regularization parameters are
  supplied

## densityratio 0.2.2

CRAN release: 2025-07-18

- update
  [`plot_bivariate()`](https://thomvolker.github.io/densityratio/reference/plot_bivariate.md)
  to depend on `ggh4x` to remove empty panels, instead of a grob (such
  that it remains a ggplot) object
- update link in readme

## densityratio 0.2.1

CRAN release: 2025-06-18

- patch test files to work with the new upcoming ggplot release
- update link readme

## densityratio 0.2.0

CRAN release: 2025-05-19

- First CRAN release of the `densityratio` package.
- Estimation methods
  [`kliep()`](https://thomvolker.github.io/densityratio/reference/kliep.md),
  [`kmm()`](https://thomvolker.github.io/densityratio/reference/kmm.md),
  [`lhss()`](https://thomvolker.github.io/densityratio/reference/lhss.md),
  [`naive()`](https://thomvolker.github.io/densityratio/reference/naive.md),
  [`spectral()`](https://thomvolker.github.io/densityratio/reference/spectral.md),
  [`ulsif()`](https://thomvolker.github.io/densityratio/reference/ulsif.md).
- Cross-validation for all methods incorporated (except naive).
- S3 methods [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) incorporated.
- Extensive checks for input data and parameters.
- Test files for all methods.
- Added a `NEWS.md` file to track changes to the package.
