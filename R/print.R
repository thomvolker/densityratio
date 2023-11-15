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
  cat("Pearson divergence between P(nu) and P(de): ", paste(format(x$PE, digits = digits, ...)), "\n", sep = "")
  if (!is.null(x$p_value)) {
    cat("Pr(P(nu)=P(de))",
        ifelse(x$p_value < 0.001,
               paste(" < .001"),
               paste(" = ", format(x$p_value, digits = 3, ...))),
        "\nBonferroni-corrected for testing with r(x) = P(nu)/P(de) AND r*(x) = P(de)/P(nu).",
        "\n\n", sep = "")
  } else {
    cat("For a two-sample homogeneity test, use 'summary(x, test = TRUE)'.\n\n")
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

#' Print a \code{summary.kliep} object
#'
#' @rdname print
#' @return \code{NULL}
#' @method print summary.kliep
#' @importFrom utils str
#' @export

print.summary.kliep <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste0(deparse(x$call)), "\n", sep = "")
  cat("\n")
  cat("Kernel Information:\n",
      "  Kernel type: Gaussian with L2 norm distances\n",
      "  Number of kernels: ", paste0(nrow(x$centers)), "\n", sep = "")

  cat("  Optimal sigma: ", paste(format(x$sigma_opt, digits, ...)), "\n",
      "  Optimal lambda: ", paste(format(x$lambda_opt, digits, ...)), "\n", sep = "")
  cat("  Optimal kernel weights (loocv):")
  cat(str(x$alpha_opt), "\n")
  cat("Kullback-Leibler divergence between P(nu) and P(de): ", paste(format(x$UKL, digits = digits, ...)), "\n", sep = "")
  if (!is.null(x$p_value)) {
    cat("Pr(P(nu)=P(de))",
        ifelse(x$p_value < 0.001,
               paste(" < .001"),
               paste(" = ", format(x$p_value, digits = 3, ...))),
        "\nBonferroni-corrected for testing with r(x) = P(nu)/P(de) AND r*(x) = P(de)/P(nu).",
        "\n\n", sep = "")
  } else {
    cat("For a two-sample homogeneity test, use 'summary(x, test = TRUE)'.\n\n")
  }
  invisible(x)
}

#' Print a \code{naivedensityratio} object
#'
#' @rdname print
#' @return \code{invisble} The inputted \code{naivedensityratio} object.
#' @method print naivedensityratio
#' @importFrom utils str
#' @export
print.naivedensityratio <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste0(deparse(x$call)), "\n", sep = "")
  cat("\n")
  cat("Naive density ratio\n",
      "  Number of variables: ", ncol(as.matrix(x$df_numerator)), "\n",
      "  Number of numerator samples: ", nrow(as.matrix(x$df_numerator)), "\n",
      "  Number of denominator samples: ", nrow(as.matrix(x$df_denominator)), "\n", sep = "")
  cat("  Numerator density:")
  cat(str(stats::predict(x, newdata = x$df_numerator)))
  cat("  Denominator density:")
  cat(str(stats::predict(x, newdata = x$df_denominator)), "\n")

  invisible(x)
}

#' Print a \code{naivesubspacedensityratio} object
#'
#' @rdname print
#' @return \code{invisble} The inputted \code{naivesubspacedensityratio} object.
#' @method print naivesubspacedensityratio
#' @importFrom utils str
#' @export
print.naivesubspacedensityratio <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste0(deparse(x$call)), "\n", sep = "")
  cat("\n")
  cat("Naive-subspace density ratio\n",
      "  Number of variables: ", ncol(as.matrix(x$df_numerator)), "\n",
      "  Size of subspace: ", x$subspace_dim, "\n",
      "  Number of numerator samples: ", nrow(as.matrix(x$df_numerator)), "\n",
      "  Number of denominator samples: ", nrow(as.matrix(x$df_denominator)), "\n", sep="")
  cat("  Numerator density:")
  cat(str(stats::predict(x, newdata = x$df_numerator)))
  cat("  Denominator density:")
  cat(str(stats::predict(x, newdata = x$df_denominator)), "\n\n")

  invisible(x)
}

#' Print a \code{summary.naivedensityratio} object
#'
#' @rdname print
#' @return \code{NULL}
#' @method print summary.naivedensityratio
#' @importFrom utils str
#' @export

print.summary.naivedensityratio <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste0(deparse(x$call)), "\n", sep = "")
  cat("\n")
  cat("Naive density ratio estimate:\n",
      "  Number of variables: ", x$nvars, "\n",
      "  Number of numerator samples: ", x$n[1], "\n",
      "  Number of denominator samples: ", x$n[2], "\n", sep="")
  cat("  Density ratio for numerator samples:")
  cat(str(x$dr$dr[1:x$n[1]]))
  cat("  Density ratio for denominator samples:")
  cat(str(x$dr$dr[(x$n[1]+1):(x$n[1]+x$n[2])]), "\n\n")

  cat("Squared average log density ratio difference for numerator and denominator samples (SALDRD): ",
      paste(format(x$SALDRD, digits = digits, ...)), "\n", sep = "")
  if (!is.null(x$p_value)) {
    cat("Pr(P(nu)=P(de))",
        ifelse(x$p_value < 0.001,
               paste(" < .001"),
               paste(" = ", format(x$p_value, digits = 3, ...))),
        "\n\n", sep = "")
  } else {
    cat("For a two-sample homogeneity test, use 'summary(x, test = TRUE)'.\n\n")
  }
  invisible(x)
}

#' Print a \code{summary.naivesubspacedensityratio} object
#'
#' @rdname print
#' @return \code{NULL}
#' @method print summary.naivesubspacedensityratio
#' @importFrom utils str
#' @export


print.summary.naivesubspacedensityratio <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste0(deparse(x$call)), "\n", sep = "")
  cat("\n")
  cat("Naive-subspace density ratio\n",
      "  Number of variables: ", x$nvars, "\n",
      "  Dimension of subspace: ", x$subspace_dim, "\n",
      "  Number of numerator samples: ", x$n[1], "\n",
      "  Number of denominiator samples: ", x$n[2], "\n", sep="")
  cat("  Density ratio for numerator samples:")
  cat(str(x$dr$dr[1:x$n[1]]))
  cat("  Density ratio for denominator samples:")
  cat(str(x$dr$dr[(x$n[1]+1):(x$n[1]+x$n[2])]), "\n\n")

  cat("Squared average log density ratio difference for numerator and denominator samples (SALDRD): ",
      paste(format(x$SALDRD, digits = digits, ...)), "\n", sep = "")
  if (!is.null(x$p_value)) {
    cat("Pr(P(nu)=P(de))",
        ifelse(x$p_value < 0.001,
               paste(" < .001"),
               paste(" = ", format(x$p_value, digits = 3, ...))),
        "\n\n", sep = "")
  } else {
    cat("For a two-sample homogeneity test, use 'summary(x, test = TRUE)'.\n\n", sep="")
  }
  invisible(x)
}
