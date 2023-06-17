#' Print a \code{ulsif} object
#'
#' @rdname print
#' @param object Object of class \code{ulsif} or \code{kliep}
#' @return \code{NULL}
#' @method print ulsif
#' @importFrom utils str
#' @export

print.ulsif <- function(object, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste0(deparse(object$call)), "\n", sep = "")
  cat("\n")
  cat("Kernel Information:\n",
      "  Kernel type: Gaussian with L2 norm distances\n",
      "  Number of kernels: ", paste0(nrow(object$centers)), "\n", sep = "")
  cat("  sigma:")
  cat(str(unname(object$sigma)))
  cat("  lambda:")
  cat(str(unname(object$lambda)))
  cat("  Optimal sigma: ", paste(format(object$sigma_opt, digits)), "\n",
      "  Optimal lambda: ", paste(format(object$lambda_opt, digits)), "\n", sep = "")
  cat("  Optimal kernel weights (loocv):")
  cat(str(object$alpha_opt), "\n")
}

#' Print a \code{ulsif} object
#'
#' @rdname print
#' @return \code{NULL}
#' @method print kliep
#' @importFrom utils str
#' @export

print.kliep <- function(object, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste0(deparse(object$call)), "\n", sep = "")
  cat("\n")
  cat("Kernel Information:\n",
      "  Kernel type: Gaussian with L2 norm distances\n",
      "  Number of kernels: ", paste0(nrow(object$centers)), "\n", sep = "")
  cat("  sigma:")
  cat(str(unname(object$sigma)))
  if (!is.null(object$cv_score)) {
    cat("  Optimal sigma (", paste(object$nfold), "-fold cv): ", paste(format(object$sigma_opt, digits = digits), collapse = "  "), "\n", sep = "")
    cat("  Optimal kernel weights (", paste(object$nfold), "-fold cv): ", sep = "")
    cat(str(object$alpha_opt))
  } else {
    cat("  Optimal sigma: NULL (no cross-validation)\n", sep = "")
    cat("  Optimal kernel weights: NULL (no cross-validation)\n", sep = "")
  }
  cat("\nOptimization parameters:\n", sep = "")
  cat("  Learning rate (epsilon): ", paste(format(object$epsilon, digits = digits), collapse = "  "), "\n", sep = "")
  cat("  Maximum number of iterations: ", paste(object$maxit, collapse = "  "))
  cat("\n")
  invisible(object)
}
