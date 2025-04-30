
test_that("distance and kernel work", {
  x <- matrix(1:10)
  y <- matrix(20:11)

  D <- distance(x, y)

  expect_equal(
    D,
    as.matrix(dist(rbind(x, y), method = "euclidean"))[11:20, 1:10] ** 2,
    ignore_attr = TRUE
  )
  expect_equal(
    kernel_gaussian(D, 2),
    exp(-D / (2*2^2)),
    ignore_attr = TRUE
  )

  expect_error(
    distance(as.data.frame(x), as.data.frame(y))
  )
  expect_error(
    kernel_gaussian(as.data.frame(D), 2)
  )
})
