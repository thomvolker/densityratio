
check.dataform <- function(nu, de) {
  if (! (is.numeric(nu) & is.numeric(de))) {
    stop("Currently only numeric data is supported.")
  }
  if (ncol(nu) != ncol(de) | !all(names(nu) == names(de))) {
    stop("nu and de must contain exactly the same set of variables.")
  }
}
