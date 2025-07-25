# densityratio 0.2.2

* update `plot_bivariate()` to depend on `ggh4x` to remove empty panels, instead
of a grob (such that it remains a ggplot) object
* update link in readme

# densityratio 0.2.1

* patch test files to work with the new upcoming ggplot release
* update link readme

# densityratio 0.2.0

* First CRAN release of the `densityratio` package.
* Estimation methods `kliep()`, `kmm()`, `lhss()`, `naive()`, `spectral()`,
`ulsif()`.
* Cross-validation for all methods incorporated (except naive).
* S3 methods `predict()`, `plot()`, `print()` and `summary()` incorporated.
* Extensive checks for input data and parameters.
* Test files for all methods.
* Added a `NEWS.md` file to track changes to the package.
