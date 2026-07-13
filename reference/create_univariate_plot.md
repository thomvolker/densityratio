# Univariate plot

Scatterplot of individual values and density ratio estimates. Used
internally in `create_univariate_plot()`

## Usage

``` r
create_univariate_plot(data, ext, var, y_lab, sample.facet = TRUE)
```

## Arguments

- data:

  Data frame with the individual values and density ratio estimates

- ext:

  Data frame with the density ratio estimates and sample indicator

- var:

  Name of the variable to be plotted on the x-axis

- y_lab:

  Name of the y-axis label, typically ("Density Ratio" or "Log Density
  Ratio")

- sample.facet:

  Logical indicating whether to facet the plot by sample. Default is
  TRUE.

## Value

A scatterplot of variable values and density ratio estimates.
