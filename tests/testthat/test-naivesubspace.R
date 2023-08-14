source(test_path("generate-test-data.R"))

test_that("1-dimensional naive subspace estimation and prediction works", {
  dr <- naivesubspace(test_df_numerator_1, test_df_denominator_1)
  expect_s3_class(dr, "naivesubspacedensityratio")
  pred <- predict(dr, test_df_newdata_1)
  expect_gt(pred[1], pred[2])
})

test_that("10-dimensional naive subspace estimation and prediction works", {
  dr <- naivesubspace(test_df_numerator_10, test_df_denominator_10)
  expect_s3_class(dr, "naivesubspacedensityratio")
  pred <- predict(dr, test_df_newdata_10)
  expect_gt(pred[1], pred[2])
})

test_that("high-dimensional naive subspace estimation and prediction works", {
  dr <- naivesubspace(test_df_numerator_max, test_df_denominator_max)
  expect_s3_class(dr, "naivesubspacedensityratio")
  pred <- predict(dr, test_df_newdata_max)
  expect_gt(pred[1], pred[2])
})
