# Bivariate plot

Bivariate plot

## Usage

``` r
create_bivariate_plot(data, ext, vars, logscale, show.sample)
```

## Arguments

- data:

  Data frame with the individual values and density ratio estimates

- ext:

  Data frame with the density ratio estimates and sample indicator

- vars:

  Character vector of variable names to be plotted.

- logscale:

  Logical indicating whether the density ratio should be plotted in log
  scale. Defaults to TRUE.

- show.sample:

  Logical indicating whether to give different shapes to observations,
  depending on the sample they come from (numerator or denominator).
  Defaults to FALSE.

## Value

Bivariate plot
