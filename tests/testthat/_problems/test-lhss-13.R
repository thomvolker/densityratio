# Extracted from test-lhss.R:13

# test -------------------------------------------------------------------------
set.seed(1)
dr <- lhss(numerator_small$x3, denominator_small$x3)
summdr <- summary(dr)
expect_s3_class(dr, "lhss")
expect_s3_class(summdr, "summary.lhss")
expect_invisible(print(summdr))
pred <- predict(dr)
expect_gt(mean(log(pmax(1e-3, pred))), 0)
expect_lt(mean(log(pmax(1e-3, predict(dr, denominator_small$x3)))), 0.01)
