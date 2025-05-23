% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{create_bivariate_plot}
\alias{create_bivariate_plot}
\title{Bivariate plot}
\usage{
create_bivariate_plot(data, ext, vars, logscale, show.sample)
}
\arguments{
\item{data}{Data frame with the individual values and density ratio estimates}

\item{ext}{Data frame with the density ratio estimates and sample indicator}

\item{vars}{Character vector of variable names to be plotted.}

\item{logscale}{Logical indicating whether the density ratio should be plotted in log scale. Defaults to TRUE.}

\item{show.sample}{Logical indicating whether to give different shapes to observations, depending on the sample they come from (numerator or denominator). Defaults to FALSE.}
}
\value{
Bivariate plot
}
\description{
Bivariate plot
}
