#' Print a \code{ulsif} object
#'
#' @rdname print
#' @param x Object of class \code{ulsif}, \code{summary.ulsif}, \code{kliep}
#' or \code{summary.kliep}.
#' @param digits Number of digits to use when printing the output.
#' @param ... further arguments on how to format the number of digits.
#' @return \code{invisble} The inputted \code{ulsif} object.
#' @method print ulsif
#' @importFrom utils str
#' @export

print.ulsif <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste0(deparse(x$call)), "\n", sep = "")
  cat("\n")
  cat("Kernel Information:\n",
      "  Kernel type: Gaussian with L2 norm distances\n",
      "  Number of kernels: ", paste0(nrow(x$centers)), "\n", sep = "")
  cat("  sigma:")
  cat(str(unname(x$sigma)))
  cat("  lambda:")
  cat(str(unname(x$lambda)))
  cat("  Optimal sigma: ", paste(format(x$sigma_opt, digits, ...)), "\n",
      "  Optimal lambda: ", paste(format(x$lambda_opt, digits, ...)), "\n", sep = "")
  cat("  Optimal kernel weights (loocv):")
  cat(str(x$alpha_opt), "\n")
  invisible(x)
}

#' Print a \code{summary.ulsif} object
#'
#' @rdname print
#' @return \code{NULL}
#' @method print summary.ulsif
#' @importFrom utils str
#' @export

print.summary.ulsif <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste0(deparse(x$call)), "\n", sep = "")
  cat("\n")
  cat("Kernel Information:\n",
      "  Kernel type: Gaussian with L2 norm distances\n",
      "  Number of kernels: ", paste0(nrow(x$centers)), "\n", sep = "")

  cat("  Optimal sigma: ", paste(format(x$sigma_opt, digits, ...)), "\n",
      "  Optimal lambda: ", paste(format(x$lambda_opt, digits, ...)), "\n", sep = "")
  cat("  Optimal kernel weights (loocv):")
  cat(str(x$alpha_opt), "\n")
  #TODO: Check pearson divergence interpretation
  cat("Pearson divergence between P(nu) and P(de): ", paste(format(x$PE, digits = digits, ...)), "\n", sep = "")
  if (!is.null(x$p_value)) {
    cat("Pr(P(nu)=P(de)) = ", paste(format(x$p_value, digits = 3, ...)), "\n\n", sep = "")
  }
  invisible(x)
}

#' Print a \code{kliep} object
#'
#' @rdname print
#' @return \code{NULL}
#' @method print kliep
#' @importFrom utils str
#' @export

print.kliep <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste0(deparse(x$call)), "\n", sep = "")
  cat("\n")
  cat("Kernel Information:\n",
      "  Kernel type: Gaussian with L2 norm distances\n",
      "  Number of kernels: ", paste0(nrow(x$centers)), "\n", sep = "")
  cat("  sigma:")
  cat(str(unname(x$sigma)))
  if (!is.null(x$cv_score)) {
    cat("  Optimal sigma (", paste(x$nfold), "-fold cv): ", paste(format(x$sigma_opt, digits = digits, ...), collapse = "  "), "\n", sep = "")
    cat("  Optimal kernel weights (", paste(x$nfold), "-fold cv): ", sep = "")
    cat(str(x$alpha_opt))
  } else {
    cat("  Optimal sigma: NULL (no cross-validation)\n", sep = "")
    cat("  Optimal kernel weights: NULL (no cross-validation)\n", sep = "")
  }
  cat("\nOptimization parameters:\n", sep = "")
  cat("  Learning rate (epsilon): ", paste(format(x$epsilon, digits = digits, ...), collapse = "  "), "\n", sep = "")
  cat("  Maximum number of iterations: ", paste(x$maxit, collapse = "  "))
  cat("\n")
  invisible(x)
}
