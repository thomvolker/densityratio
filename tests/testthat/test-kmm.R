
test_that("1-dimensional kmm estimation and prediction works", {
  set.seed(1)
  dr <- kmm(numerator_small$x3, denominator_small$x3)
  summdr <- summary(dr)
  expect_s3_class(dr, "kmm")
  expect_s3_class(summdr, "summary.kmm")
  expect_invisible(print(summdr))

  expect_equal(summdr$centers, dr$centers)
  expect_equal(summdr$alpha_opt, dr$alpha_opt)
  expect_equal(summdr$sigma_opt, dr$sigma_opt)

  pred <- predict(dr)[, 1]
  expect_gt(mean(log(pmax(1e-3, pred))), 0)
  expect_lt(mean(log(pmax(1e-3, predict(dr, denominator_small$x3)[,1]))), 0)

  dr <- kmm(numerator_small$x3, denominator_small$x3, sigma = 2, scale = NULL)
  summdr <- summary(dr, test = TRUE)
  expect_lte(summdr$p_value, 1)

  expect_invisible(print(dr))
  expect_invisible(print(summdr))
})

test_that("multidimensional kmm estimation, prediction and plotting works", {
  set.seed(1)
  dr <- kmm(numerator_small, denominator_small, progressbar = FALSE)
  expect_s3_class(dr, "kmm")

  expect_true(
    ggplot2::ggplot2::is_ggplot(
      suppressWarnings(plot(dr))
    )
  )

  expect_gt(mean(log(pmax(1e-3, predict(dr)))), 0)
  expect_lt(mean(log(pmax(1e-3, predict(dr, denominator_small)[,1]))), 0)

  dr <- kmm(numerator_small, denominator_small, sigma = c(2, 3), ncenters = 100, scale = NULL,
            constrained = TRUE)

  expect_type(dr$alpha_opt, "double")
  expect_type(dr$sigma, "double")

  Kdn <- distance(
    as.matrix(denominator_small),
    as.matrix(numerator_small),
    FALSE
  ) |> kernel_gaussian(2)

  Kdd <- distance(
    as.matrix(denominator_small),
    as.matrix(denominator_small),
    FALSE
  ) |> kernel_gaussian(2)

  est <- compute_kmm(
    as.matrix(numerator_small),
    as.matrix(denominator_small),
    as.matrix(denominator_small),
    distance(as.matrix(denominator_small), as.matrix(denominator_small)),
    2,
    rep(0, nrow(numerator_small)),
    rep(0, nrow(denominator_small)),
    FALSE,
    1,
    FALSE,
    FALSE,
    NULL
  )

  expect_equal(
    est$alpha,
    solve(crossprod(Kdd) %*% Kdd + 1e-3*diag(nrow(Kdd)),
          rowSums(crossprod(Kdd, Kdn))) *
      nrow(Kdd) / ncol(Kdn),
    ignore_attr = TRUE
  )
  expect_equal(est$loss, 0, ignore_attr = TRUE)

  est_constrained <- compute_kmm(
    as.matrix(numerator_small),
    as.matrix(denominator_small),
    as.matrix(denominator_small),
    distance(as.matrix(denominator_small), as.matrix(denominator_small)),
    2,
    rep(0, ncol(numerator_small)),
    rep(0, nrow(denominator_small)),
    FALSE,
    1,
    FALSE,
    TRUE,
    NULL
  )

  alpha_constrained <- kmm_constrained_alpha(
    Kdn, Kdd, Kdd, ncol(Kdn), nrow(Kdn), NULL
  )

  expect_gte(min(est_constrained$alpha) |> round(2), 0)
  expect_equal(est$loss, 0, ignore_attr = TRUE)
  expect_equal(
    est_constrained$alpha,
    alpha_constrained,
    ignore_attr = TRUE
  )
})
