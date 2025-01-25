#' numerator_data
#'
#' Simulated data set (see data-raw/generate-data-densityratio.R) with five
#' variables that are used in the examples.
#'
#' @name numerator_data
#' @docType data
#' @format A data frame with 1000 rows and 5 columns:
#' \describe{
#'    \item{x1}{Categorical variable with three categories, 'A', 'B' and 'C'}
#'    \item{x2}{Categorical variable with two categories, 'G1' and 'G2'}
#'    \item{x3}{Continuous variable (normally distributed given x1 and x2)}
#'    \item{x4}{Continuous variable (normally distributed given x3)}
#'    \item{x5}{Continuous variable (mixture of two normally distributed variables)}
#' }
#' @keywords data datasets

NULL

#' denominator_data
#'
#' Simulated data set (see data-raw/generate-data-densityratio.R) with five
#' variables that are used in the examples.
#'
#' @name denominator_data
#' @docType data
#' @format A data frame with 1000 rows and 5 columns:
#' \describe{
#'    \item{x1}{Categorical variable with three categories, 'A', 'B' and 'C'}
#'    \item{x2}{Categorical variable with two categories, 'G1' and 'G2'}
#'    \item{x3}{Continuous variable (normally distributed given x1 and x2)}
#'    \item{x4}{Continuous variable (normally distributed)}
#'    \item{x5}{Continuous variable (normally distributed)}
#' }
#' @keywords data datasets

NULL

#' numerator_small
#'
#' Subset of the [numerator_data] with three variables and 50 observations
#'
#' @name numerator_small
#' @docType data
#' @format A data frame with 50 rows and 3 columns:
#' \describe{
#'    \item{x1}{Continuous variable (normally distributed given x1 and x2)}
#'    \item{x2}{Continuous variable (normally distributed given x3)}
#'    \item{x3}{Continuous variable (mixture of two normally distributed variables)}
#' }
#' @keywords data datasets

NULL

#' denominator_small
#'
#' Subset of the [denominator_data] with three variables and 50 observations
#'
#' @name denominator_small
#' @docType data
#' @format A data frame with 100 rows and 3 columns:
#' \describe{
#'    \item{x1}{Continuous variable (normally distributed given x1 and x2)}
#'    \item{x2}{Continuous variable (normally distributed)}
#'    \item{x3}{Continuous variable (normally distributed)}
#' }
#' @keywords data datasets

NULL



