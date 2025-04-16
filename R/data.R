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

#' insurance
#'
#' Insurance data that is openly available (e.g., on
#' \href{Kaggle}{https://www.kaggle.com/datasets/mirichoi0218/insurance}).
#'
#' @name insurance
#' @docType data
#' @format A data.frame with 1338 rows and 7 columns:
#' \describe{
#'   \item{age}{Age of the insured (continuous)}
#'   \item{sex}{Sex of the insured (binary)}
#'   \item{bmi}{Body mass index of the insured (continuous)}
#'   \item{children}{Number of children/dependents covered by the insurance (integer)}
#'   \item{smoker}{Whether the insured is a smoker (binary)}
#'   \item{region}{The region in which the insured lives (categorical)}
#'   \item{charges}{The medical costs billed by the insurance (continuous)}
#' }
#' @keywords data datasets insurance

NULL

#' colon
#'
#' Colon cancer data set from princeton, containing 2000 gene expressions from
#' 22 colon tumor tissues and 40 non-tumor tissues. The data can also be found
#' \href{here}{http://genomics-pubs.princeton.edu/oncology/affydata/I2000.html}.
#'
#' @name colon
#' @docType data
#' @format A data.frame with 62 rows and 2001 columns (class variable and 2000
#' gene expressions).
#' @keywords  data datasets colon
#' @references
#' Alon, U., Barkai, N., Notterman, D. A., Gish, K., Ybarra, S., Mack, D., &
#' Levine, A. J. (1999). Broad patterns of gene expression revealed by clustering
#' of tumor and normal colon tissues probed by oligonucleotide arrays.
#' *Proceedings of the National Academy of Sciences*, 96(12), 6745-6750.

NULL

#' kidiq
#'
#' The kidiq data stems from the National Longitudinal Survey of Youth and is used
#' in Gelman and Hill (2007). The data set contains 434 observations measured on
#' five variables, and is obtained from
#' \href{https://github.com/jknowles/BDAexampleR}{https://github.com/jknowles/BDAexampleR}.
#'
#' @name kidiq
#' @docType data
#' @format A data.frame with 434 rows and 5 columns
#' \describe{
#'   \item{kid_score}{Child's IQ score (continuous)}
#'   \item{mom_hs}{Whether the mother obtained a high school degree (binary)}
#'   \item{mom_iq}{Mother's IQ score (continuous)}
#'   \item{mom_work}{Whether the mother worked in the first three years of the
#'   child's life (1: not in the first three years; 2: in the second or third
#'   year; 3: parttime in the first year; 4: fulltime in the first year)}
#'   \item{mom_age}{Mother's age (continuous)}
#' }
#' @keywords data datasets kidiq
#' @references
#' Gelman, A., & Hill, J. (2006). Data Analysis Using Regression and
#' Multilevel/Hierarchical Models. Cambridge: Cambridge University Press.

NULL
