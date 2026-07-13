# Extracted from test-lhss.R:37

# test -------------------------------------------------------------------------
set.seed(1)
dr <- lhss(numerator_small, denominator_small, m = 2, sigma = c(0.1, 1, 2, 3),
             lambda = c(1, 0.1, 0.001),
             intercept = FALSE)
expect_s3_class(dr, "lhss")
expect_true(
    ggplot2::is_ggplot(
      suppressWarnings(plot(dr))
    )
  )
expect_gt(mean(log(pmax(1e-3, predict(dr)))), 0)
expect_lt(mean(log(pmax(1e-3, predict(dr, denominator_small)))), 0)
