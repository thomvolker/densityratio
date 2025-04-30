
test_that("1-dimensional kliep estimation and prediction works", {
  set.seed(1)
  dr <- kliep(numerator_small$x3, denominator_small$x3)
  summdr <- summary(dr)
  expect_s3_class(dr, "kliep")
  expect_s3_class(summdr, "summary.kliep")
  expect_invisible(print(summdr))

  expect_equal(summdr$centers, dr$centers)
  expect_equal(summdr$alpha_opt, dr$alpha_opt)
  expect_equal(summdr$sigma_opt, dr$sigma_opt)

  pred <- predict(dr)[, 1]
  expect_gt(mean(log(pmax(1e-3, pred))), 0)
  expect_lt(mean(log(pmax(1e-3, predict(dr, denominator_small$x3)[,1]))), 0)

  dr <- kliep(numerator_small$x3, denominator_small$x3,
              sigma = 2, scale = NULL)
  summdr <- summary(dr, test = TRUE)
  expect_lte(summdr$p_value, 1)

  expect_invisible(print(dr))
  expect_invisible(print(summdr))
})

test_that("multidimensional kliep estimation and prediction works", {
  set.seed(1)
  dr <- kliep(numerator_small, denominator_small)
  expect_s3_class(dr, "kliep")
  expect_type(plot(dr) |> suppressWarnings(), "list")

  expect_gt(mean(log(pmax(1e-3, predict(dr)))), 0)
  expect_lt(mean(log(pmax(1e-3, predict(dr, denominator_small)[,1]))), 0)

  dr <- kliep(numerator_small, denominator_small, sigma = 2,
              ncenters = 100, scale = NULL)

  expect_type(dr$alpha_opt, "double")
  expect_type(dr$sigma, "double")
  expect_type(predict(dr, sigma = 3), "double")

  Dnu <- distance(
    as.matrix(numerator_small),
    as.matrix(numerator_small),
    FALSE
  )
  Dde <- distance(
    as.matrix(denominator_small),
    as.matrix(numerator_small),
    FALSE
  )

  est <- compute_kliep(Dnu, Dde, 2, check.epsilon(NULL), 5000, rep(0, nrow(Dnu)), FALSE)
  expect_equal(est$cv_score, 0, ignore_attr = TRUE)
  expect_type(est$alpha, "double")
})
