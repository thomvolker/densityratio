#' Print a \code{ulsif} object
#'
#' @rdname print.ulsif
#' @method print ulsif
#' @param x Object of class \code{ulsif}.
#' @param digits Number of digits to use when printing the output.
#' @param ... further arguments on how to format the number of digits.
#' @return \code{invisble} The inputted \code{ulsif} object.
#' @importFrom utils str
#' @export
#' @seealso \code{\link{print}}, \code{\link{ulsif}}


print.ulsif <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste0(deparse(x$call)), "\n", sep = "")
  cat("\n")
  cat("Kernel Information:\n",
      "  Kernel type: Gaussian with L2 norm distances\n",
      "  Number of kernels: ", paste0(nrow(x$centers)), "\n", sep = "")
  cat("  sigma:")
  cat(str(unname(x$sigma)))
  cat("\nRegularization parameter (lambda):")
  cat(str(unname(x$lambda)))
  cat("\nOptimal sigma: ", paste(format(x$sigma_opt, digits, ...)), "\n",
      "Optimal lambda: ", paste(format(x$lambda_opt, digits, ...)), "\n", sep = "")
  cat("Optimal kernel weights (loocv):")
  cat(str(x$alpha_opt), "\n")
  invisible(x)
}

#' Print a \code{summary.ulsif} object
#'
#' @rdname print.summary.ulsif
#' @method print summary.ulsif
#' @param x Object of class \code{summary.ulsif}.
#' @param digits Number of digits to use when printing the output.
#' @param ... further arguments on how to format the number of digits.
#' @return \code{invisble} The inputted \code{summary.ulsif} object.
#' @importFrom utils str
#' @seealso \code{\link{print}}, \code{\link{summary.ulsif}}, \code{\link{ulsif}}
#'
#' @export

print.summary.ulsif <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste0(deparse(x$call)), "\n", sep = "")
  cat("\n")
  cat("Kernel Information:\n",
      "  Kernel type: Gaussian with L2 norm distances\n",
      "  Number of kernels: ", paste0(nrow(x$centers)), "\n", sep = "")

  cat("\nOptimal sigma: ", paste(format(x$sigma_opt, digits, ...)), "\n",
      "Optimal lambda: ", paste(format(x$lambda_opt, digits, ...)), "\n", sep = "")
  cat("Optimal kernel weights (loocv):")
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
#' @rdname print.kliep
#' @method print kliep
#' @param x Object of class \code{kliep}.
#' @param digits Number of digits to use when printing the output.
#' @param ... further arguments on how to format the number of digits.
#' @return \code{invisble} The inputted \code{kliep} object.
#' @importFrom utils str
#' @export
#' @seealso \code{\link{print}}, \code{\link{kliep}}

print.kliep <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste0(deparse(x$call)), "\n", sep = "")
  cat("\n")
  cat("Kernel Information:\n",
      "  Kernel type: Gaussian with L2 norm distances\n",
      "  Number of kernels: ", paste0(nrow(x$centers)), "\n", sep = "")
  cat("  sigma:")
  cat(str(unname(x$sigma)))
  if (!is.null(x$cv_score)) {
    cat("\nOptimal sigma (", paste(x$nfold), "-fold cv): ", paste(format(x$sigma_opt, digits = digits, ...), collapse = "  "), "\n", sep = "")
    cat("Optimal kernel weights (", paste(x$nfold), "-fold cv): ", sep = "")
    cat(str(x$alpha_opt))
  } else {
    cat("\nOptimal sigma: NULL (no cross-validation)\n", sep = "")
    cat("Optimal kernel weights: NULL (no cross-validation)\n", sep = "")
  }
  cat("\nOptimization parameters:\n", sep = "")
  cat("  Learning rate (epsilon): ", paste(format(x$epsilon, digits = digits, ...), collapse = "  "), "\n", sep = "")
  cat("  Maximum number of iterations: ", paste(x$maxit, collapse = "  "))
  cat("\n")
  invisible(x)
}

#' Print a \code{summary.kliep} object
#'
#' @rdname print.summary.kliep
#' @method print summary.kliep
#' @param x Object of class \code{summary.kliep}.
#' @param digits Number of digits to use when printing the output.
#' @param ... further arguments on how to format the number of digits.
#' @return \code{invisble} The inputted \code{summary.kliep} object.
#' @importFrom utils str
#' @seealso \code{\link{print}}, \code{\link{summary.kliep}}, \code{\link{kliep}}
#'
#' @export

print.summary.kliep <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste0(deparse(x$call)), "\n", sep = "")
  cat("\n")
  cat("Kernel Information:\n",
      "  Kernel type: Gaussian with L2 norm distances\n",
      "  Number of kernels: ", paste0(nrow(x$centers)), "\n", sep = "")

  cat("Optimal sigma: ", paste(format(x$sigma_opt, digits, ...)), "\n",
      "Optimal lambda: ", paste(format(x$lambda_opt, digits, ...)), "\n", sep = "")
  cat("Optimal kernel weights (loocv):")
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

#' Print a \code{lhss} object
#'
#' @rdname print.lhss
#' @method print lhss
#' @param x Object of class \code{lhss}.
#' @param digits Number of digits to use when printing the output.
#' @param ... further arguments on how to format the number of digits.
#' @return \code{invisble} The inputted \code{lhss} object.
#' @importFrom utils str
#' @export
#' @seealso \code{\link{print}}, \code{\link{lhss}}

print.lhss <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste0(deparse(x$call)), "\n", sep = "")
  cat("\n")
  cat("Kernel Information:\n",
      "  Kernel type: Gaussian with L2 norm distances\n",
      "  Number of kernels: ", paste0(nrow(x$centers)), "\n", sep = "")
  cat("  sigma:")
  cat(str(unname(x$sigma)))
  cat("\nRegularization parameter (lambda):")
  cat(str(unname(x$lambda)))
  cat("\nSubspace dimension (m): ", x$m, "\n", sep = "")
  cat("Optimal sigma: ", paste(format(x$sigma_opt, digits, ...)), "\n",
      "Optimal lambda: ", paste(format(x$lambda_opt, digits, ...)), "\n", sep = "")
  cat("Optimal kernel weights (loocv):")
  cat(str(x$alpha_opt), "\n")
  invisible(x)
}

#' Print a \code{summary.lhss} object
#'
#' @rdname print.summary.lhss
#' @method print summary.lhss
#' @param x Object of class \code{summary.lhss}.
#' @param digits Number of digits to use when printing the output.
#' @param ... further arguments on how to format the number of digits.
#' @return \code{invisble} The inputted \code{summary.lhss} object.
#' @importFrom utils str
#' @seealso \code{\link{print}}, \code{\link{summary.lhss}}, \code{\link{lhss}}
#'
#' @export

print.summary.lhss <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste0(deparse(x$call)), "\n", sep = "")
  cat("\n")
  cat("Kernel Information:\n",
      "  Kernel type: Gaussian with L2 norm distances\n",
      "  Number of kernels: ", paste0(nrow(x$centers)), "\n", sep = "")
  cat("\nSubspace dimension (m): ", x$m, "\n", sep = "")
  cat("Optimal sigma: ", paste(format(x$sigma_opt, digits, ...)), "\n",
      "Optimal lambda: ", paste(format(x$lambda_opt, digits, ...)), "\n", sep = "")
  cat("Optimal kernel weights (loocv):")
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

#' Print a \code{spectral} object
#'
#' @rdname print.spectral
#' @method print spectral
#' @param x Object of class \code{spectral}.
#' @param digits Number of digits to use when printing the output.
#' @param ... further arguments on how to format the number of digits.
#' @return \code{invisble} The inputted \code{spectral} object.
#' @importFrom utils str
#' @export
#' @seealso \code{\link{print}}, \code{\link{spectral}}


print.spectral <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste0(deparse(x$call)), "\n", sep = "")
  cat("\n")
  cat("Kernel Information:\n",
      "  Kernel type: Gaussian with L2 norm distances\n",
      "  Number of kernels: ", paste0(nrow(x$centers)), "\n", sep = "")
  cat("  sigma:")
  cat(str(unname(x$sigma)))
  cat("\n")
  cat("\nSubspace dimension (J):")
  cat(str(unname(x$J)))
  cat("\nOptimal sigma: ", paste(format(x$sigma_opt, digits, ...)), "\n",
      "Optimal subspace: ", paste(format(x$J_opt, digits, ...)), "\n", sep = "")
  cat("Optimal kernel weights (loocv):")
  cat(str(x$alpha_opt), "\n")
  invisible(x)
}

#' Print a \code{summary.spectral} object
#'
#' @rdname print.summary.spectral
#' @method print summary.spectral
#' @param x Object of class \code{summary.spectral}.
#' @param digits Number of digits to use when printing the output.
#' @param ... further arguments on how to format the number of digits.
#' @return \code{invisble} The inputted \code{summary.spectral} object.
#' @importFrom utils str
#' @seealso \code{\link{print}}, \code{\link{summary.spectral}}, \code{\link{spectral}}
#'
#' @export

print.summary.spectral <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste0(deparse(x$call)), "\n", sep = "")
  cat("\n")
  cat("Kernel Information:\n",
      "  Kernel type: Gaussian with L2 norm distances\n",
      "  Number of kernels: ", paste0(nrow(x$centers)), "\n", sep = "")

  cat("\nOptimal sigma: ", paste(format(x$sigma_opt, digits, ...)), "\n",
      "Optimal subspace: ", paste(format(x$J_opt, digits, ...)), "\n", sep = "")
  cat("Optimal kernel weights (loocv):")
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

#' Print a \code{naivedensityratio} object
#'
#' @rdname print.naivedensityratio
#' @method print naivedensityratio
#' @param x Object of class \code{naivesubspacedensityratio}.
#' @param digits Number of digits to use when printing the output.
#' @param ... further arguments on how to format the number of digits.
#' @return \code{invisble} The inputted \code{naivedensityratio} object.
#' @importFrom utils str
#' @export
#' @seealso \code{\link{print}}, \code{\link{naive}}


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


#' Print a \code{naivedensityratio} object
#'
#' @rdname print.naivesubspacedensityratio
#' @method print naivesubspacedensityratio
#' @param x Object of class \code{naivesubspacedensityratio}.
#' @param digits Number of digits to use when printing the output.
#' @param ... further arguments on how to format the number of digits.
#' @return \code{invisble} The inputted \code{naivesubspacedensityratio} object.
#' @importFrom utils str
#' @export
#' @seealso \code{\link{print}}, \code{\link{naivesubspace}}


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
#' @rdname print.summary.naivedensityratio
#' @method print summary.naivedensityratio
#' @param x Object of class \code{summary.naivedensityratio}.
#' @param digits Number of digits to use when printing the output.
#' @param ... further arguments on how to format the number of digits.
#' @return \code{invisble} The inputted \code{summary.naivedensityratio} object.
#' @importFrom utils str
#' @seealso \code{\link{print}}, \code{\link{summary.naivedensityratio}},
#' \code{\link{naive}}


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
#' @rdname print.summary.naivesubspacedensityratio
#' @method print summary.naivesubspacedensityratio
#' @param x Object of class \code{summary.naivesubspacedensityratio}.
#' @param digits Number of digits to use when printing the output.
#' @param ... further arguments on how to format the number of digits.
#' @return \code{invisble} The inputted \code{summary.naivesubspacedensityratio} object.
#' @importFrom utils str
#' @seealso \code{\link{print}}, \code{\link{summary.naivesubspacedensityratio}},
#' \code{\link{naive}}


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
