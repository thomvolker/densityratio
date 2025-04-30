
test_that("1-dimensional spectral estimation and prediction works", {
  set.seed(1)
  dr <- spectral(numerator_small$x3, denominator_small$x3)
  summdr <- summary(dr)
  expect_s3_class(dr, "spectral")
  expect_s3_class(summdr, "summary.spectral")
  expect_invisible(print(summdr))

  expect_equal(summdr$centers, dr$centers)
  expect_equal(summdr$alpha_opt, dr$alpha_opt)
  expect_equal(summdr$sigma_opt, dr$sigma_opt)

  pred <- predict(dr)[, , 1]
  expect_gt(mean(log(pmax(1e-3, pred))), 0)
  expect_lt(mean(log(pmax(1e-3, predict(dr, denominator_small$x3)[,,1]))), 0)

  dr <- spectral(numerator_small$x3, denominator_small$x3,
                 sigma = 2, scale = NULL)
  summdr <- summary(dr, test = TRUE)
  expect_lte(summdr$p_value, 1)

  expect_invisible(print(dr))
  expect_invisible(print(summdr))
})

test_that("multidimensional spectral estimation and prediction works", {
  set.seed(1)
  dr <- spectral(numerator_small, denominator_small)
  expect_s3_class(dr, "spectral")

  expect_gt(mean(log(pmax(1e-3, predict(dr)))), 0)
  expect_lt(mean(log(pmax(1e-3, predict(dr, denominator_small)[,,1]))), 0)

  dr <- spectral(numerator_small, denominator_small, m = 10, sigma = 2,
                 ncenters = 100, scale = NULL)

  expect_type(dr$alpha_opt, "double")
  expect_type(dr$sigma, "double")
  expect_type(dr$m_opt, "double")

  expect_type(
    predict(dr, sigma = pi * c(1,2,3), progressbar = FALSE),
    "double"
  )
  expect_type(
    predict(dr, sigma = pi, m = 20, progressbar = FALSE),
    "double"
  )

  Dnu <- distance(
    as.matrix(numerator_small),
    as.matrix(denominator_small),
    FALSE
  )
  Dde <- distance(
    as.matrix(denominator_small),
    as.matrix(denominator_small),
    FALSE
  )

  est <- spectral_dre(Dnu, Dde, m = 10, sigma = c(1,2,3),
                      cv_ind_nu = check.nfold(TRUE, 10, nrow(Dnu)),
                      cv_ind_de = check.nfold(TRUE, 10, nrow(Dde)),
                      parallel = FALSE, nthreads = 0, progressbar = FALSE)
  expect_type(est$betatilde, "double")
  expect_type(est$loss, "double")
  expect_type(est$Evals, "double")
  expect_type(est$Evecs, "double")
})
