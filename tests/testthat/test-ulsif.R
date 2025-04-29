
test_that("1-dimensional ULSIF estimation and prediction works", {
  set.seed(1)
  dr <- ulsif(numerator_small$x3, denominator_small$x3)
  summdr <- summary(dr)
  expect_s3_class(dr, "ulsif")
  expect_s3_class(summdr, "summary.ulsif")
  expect_invisible(print(summdr))

  expect_equal(summdr$centers, dr$centers)
  expect_equal(summdr$alpha_opt, dr$alpha_opt)
  expect_equal(summdr$sigma_opt, dr$sigma_opt)
  expect_equal(summdr$lambda_opt, dr$lambda_opt)

  pred <- predict(dr)[, , 1]
  expect_gt(mean(log(pmax(1e-3, pred))), 0)
  expect_lt(mean(log(pmax(1e-3, predict(dr, denominator_small$x3)[,,1]))), 0)

  dr <- ulsif(numerator_small$x3, denominator_small$x3, intercept = FALSE,
              lambda = 0.1, sigma = 2, centers = numerator_small$x3,
              scale = NULL)
  summdr <- summary(dr, test = TRUE)
  summdr_parallel <- summary(dr, test = TRUE, parallel = TRUE, cluster = 2)
  expect_lte(summdr$p_value, 1)
  expect_type(summdr_parallel, "list")
  expect_lte(summdr_parallel$p_value, 1)

  expect_invisible(print(dr))
  expect_invisible(print(summdr))
  expect_invisible(print(summdr_parallel))

  Knu <- distance(
    as.matrix(numerator_small$x3),
    as.matrix(numerator_small$x3)
  ) |> kernel_gaussian(2)
  Kde <- distance(
    as.matrix(denominator_small$x3),
    as.matrix(numerator_small$x3)
  ) |> kernel_gaussian(2)

  expect_equal(solve(crossprod(Kde)/nrow(Kde) + 0.1 * diag(ncol(Kde)), colMeans(Knu)), dr$alpha_opt)

})

test_that("multidimensional ULSIF estimation, prediction works", {
  set.seed(1)
  dr <- ulsif(numerator_small, denominator_small)
  expect_s3_class(dr, "ulsif")

  expect_gt(mean(predict(dr)[,,1]), 1)
  expect_lt(mean(predict(dr, denominator_small)[,,1]), 1)

  dr <- ulsif(numerator_small, denominator_small, intercept = FALSE,
              lambda = 0.1, sigma = 2, centers = numerator_small,
              scale = NULL)

  expect_type(
    plot(dr),
    "list"
  ) |> suppressWarnings()
  expect_type(
    plot(dr, samples = "numerator"),
    "list"
  ) |> suppressWarnings()
  expect_type(
    plot(dr, binwidth = 0.5),
    "list"
  ) |> suppressWarnings()
  expect_type(
    plot(dr, bins = 30),
    "list"
  ) |> suppressWarnings()
  expect_no_warning(
    plot(dr, logscale = FALSE)
  )

  expect_type(
    plot_univariate(dr),
    "list"
  ) |> suppressWarnings()

  expect_type(
    plot_univariate(dr, vars = c("x1", "x2"), samples = "denominator",
                    logscale = FALSE, grid = TRUE),
    "list"
  )
  expect_type(
    plot_univariate(dr, grid = TRUE, sample.facet = TRUE,
                    logscale = FALSE, nrow.panel = 3),
    "list"
  )
  expect_error(
    plot_univariate(dr, vars = c("a", "b"))
  )
  expect_type(
    plot_bivariate(dr),
    "list"
  ) |> suppressWarnings()
  expect_type(
    plot_bivariate(dr, vars = c("x1", "x2", "x3"), samples = "numerator",
                   logscale = FALSE, grid = TRUE),
    "list"
  )
  expect_error(
    plot_bivariate(dr, vars = c("x1", "b"))
  )

  expect_type(predict(dr, sigma = 3), "double")

  Knu <- distance(
    model.matrix(~., numerator_small),
    model.matrix(~., numerator_small)
  ) |> kernel_gaussian(2)
  Kde <- distance(
    model.matrix(~., denominator_small),
    model.matrix(~., numerator_small)
  ) |> kernel_gaussian(2)

  expect_equal(solve(crossprod(Kde)/nrow(Kde) + 0.1 * diag(ncol(Kde)), colMeans(Knu)), dr$alpha_opt)
})

test_that("ULSIF estimation functions work", {
  Dnu <- distance(
    model.matrix(~., numerator_small),
    model.matrix(~., numerator_small)
  )
  Dde <- distance(
    model.matrix(~., denominator_small),
    model.matrix(~., numerator_small)
  )

  ulsif_out <- compute_ulsif(
    Dnu,
    Dde,
    sigma = 2,
    lambda = 0.1,
    parallel = FALSE,
    nthreads = 1,
    progressbar = FALSE
  )

  expect_length(ulsif_out, 2)
  expect_type(ulsif_out$alpha, "double")
  expect_type(ulsif_out$loocv_score, "double")
})

test_that("set_threads works", {
  expect_equal(set_threads(10), 10)
  expect_equal(set_threads(0), parallel::detectCores())
  expect_warning(set_threads(50))
})
