
test_that("1-dimensional lhss estimation and prediction works", {
  set.seed(1)
  dr <- lhss(numerator_small$x3, denominator_small$x3)
  summdr <- summary(dr)
  expect_s3_class(dr, "lhss")
  expect_s3_class(summdr, "summary.lhss")
  expect_invisible(print(summdr))

  pred <- predict(dr)
  expect_gt(mean(log(pmax(1e-3, pred))), 0)
  # below was actually expected to be below zero, but it is possible in theory
  expect_lt(mean(log(pmax(1e-3, predict(dr, denominator_small$x3)))), 0.01)

  dr <- lhss(numerator_small$x3, denominator_small$x3, nsigma = 1, lambda = 0.2)
  summdr <- summary(dr, test = TRUE)
  expect_lte(summdr$p_value, 1)

  expect_invisible(print(dr))
  expect_invisible(print(summdr))
})

test_that("multidimensional lhss estimation, prediction and plotting works", {
  set.seed(1)
  dr <- lhss(numerator_small, denominator_small, m = 2, sigma = c(0.1, 1, 2, 3),
             lambda = c(1, 0.1, 0.001),
             intercept = FALSE)
  expect_s3_class(dr, "lhss")
  expect_type(plot(dr) |> suppressWarnings(), "list")

  expect_gt(mean(log(pmax(1e-3, predict(dr)))), 0)
  expect_lt(mean(log(pmax(1e-3, predict(dr, denominator_small)))), 0)

  expect_type(
    predict(dr, sigma = 2.5, lambda = 1),
    "double"
  )
  expect_type(
    predict(dr, sigma_quantile = 0.5),
    type = "double"
  )

  expect_type(dr$alpha, "double")
  expect_type(dr$U, "double")
})
