
check.dataform <- function(nu, de) {
  if (! (is.numeric(nu) & is.numeric(de))) {
    stop("Currently only numeric data is supported.")
  }
  if (ncol(nu) != ncol(de) | !all(names(nu) == names(de))) {
    stop("nu and de must contain exactly the same set of variables.")
  }
}

check.sigma <- function(sigma) {
  if (!(is.null(sigma) | is.numeric(sigma)) | !is.null(dim(sigma))) {
    stop("sigma must be either NULL, a numeric scalar or a numeric vector")
  }
}


check.lambda <- function(lambda) {
  if (!(is.null(lambda) | is.numeric(lambda)) | !is.null(dim(lambda))) {
    stop("lambda must be either NULL, a numeric scalar or a numeric vector")
  }
}
