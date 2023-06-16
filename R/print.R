#' @importFrom utils str
print.ulsif <- function(object) {
  cat("\nCall:\n", paste0(deparse(object$call)), "\n", sep = "")
  cat("\n")
  cat("Kernel Information:\n",
      "  Kernel type: Gaussian with L2 norm distances\n",
      "  Number of kernels: ", paste0(nrow(object$centers)), "\n", sep = "")
  cat("  sigma:")
  cat(str(unname(object$sigma)))
  cat("  lambda:")
  cat(str(unname(object$lambda)))
  cat("  Optimal sigma: ", paste0(signif(object$sigma_opt, 3)), "\n",
      "  Optimal lambda: ", paste0(signif(object$lambda_opt, 3)), "\n", sep = "")
  cat("  Optimal kernel weights (loocv):")
  cat(str(object$alpha_opt), "\n")
}


print.kliep <- function(object) {
  cat("\nCall:\n", paste0(deparse(object$call)), "\n", sep = "")
  cat("\n")
  cat("Kernel Information:\n",
      "  Kernel type: Gaussian with L2 norm distances\n",
      "  Number of kernels: ", paste0(nrow(object$centers)), "\n", sep = "")
  cat("  sigma:")
  cat(str(unname(object$sigma)))
  if (!is.null(object$cv_score)) {
    cat("  Optimal sigma (", paste0(object$nfold), "-fold cv): ", paste0(signif(object$sigma_opt, 3)), "\n", sep = "")
    cat("  Optimal kernel weights (", paste0(object$nfold), "-fold cv): ", sep = "")
    cat(str(object$alpha_opt))
  } else {
    cat("  Optimal sigma: NULL (no cross-validation)\n", sep = "")
    cat("  Optimal kernel weights: NULL (no cross-validation)\n", sep = "")
  }
  cat("\nOptimization parameters:\n", sep = "")
  cat("  Learning rate (epsilon): ", paste(object$epsilon, collapse = "  "), "\n", sep = "")
  cat("  Maximum number of iterations: ", paste(object$maxit, collapse = "  "))
}
