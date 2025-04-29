
test_that("1-dimensional naive estimation and prediction works", {
  dr <- naive(numerator_small$x3, denominator_small$x3)
  summdr <- summary(dr)
  expect_s3_class(dr, "naivedensityratio")
  expect_s3_class(summdr, "summary.naivedensityratio")
  expect_invisible(print(summdr))

  pred <- predict(dr)
  expect_gt(mean(log(pmax(1e-3, pred))), 0)
  expect_lt(mean(log(pmax(1e-3, predict(dr, denominator_small$x3)))), 0)

  dr <- naive(numerator_small$x3, denominator_small$x3, kernel = "ep")
  summdr <- summary(dr, test = TRUE)
  expect_lte(summdr$p_value, 1)

  expect_invisible(print(dr))
  expect_invisible(print(summdr))
})

test_that("multidimensional naive estimation and prediction works", {
  dr <- naive(numerator_small, denominator_small)
  expect_s3_class(dr, "naivedensityratio")
  expect_type(plot(dr), "list")

  expect_gt(mean(log(pmax(1e-3, predict(dr)))), 0)
  expect_lt(mean(log(pmax(1e-3, predict(dr, denominator_small)))), 0)

  dr <- naive(numerator_small, denominator_small, m = 2, kernel = "cos", bw = "nrd0")

  expect_type(dr$density_numerator$PC1$y, "double")
})
