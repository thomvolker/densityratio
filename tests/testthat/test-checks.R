
test_that("check.datatype works", {
  expect_s3_class(
    check.datatype(1:10),
    "data.frame"
  )
  expect_s3_class(
    check.datatype(data.frame(x = 1:10, y = 1:10)),
    "data.frame"
  )
  expect_s3_class(
    check.datatype(matrix(1:10, 5)),
    "data.frame"
  )
  expect_error(
    check.datatype(c(1, 2, NA, 4))
  )
})

test_that("check.dataform works", {
  expect_type(
    check.dataform(
      data.frame(x = 1:10),
      data.frame(x = 10:1),
      NULL,
      TRUE,
      scale = NULL
    ),
    "list"
  )
  expect_error(
    check.dataform(
      data.frame(x = 1:10),
      data.frame(y = 10:1),
      NULL,
      TRUE,
      scale = NULL
    )
  )
  expect_error(
    check.dataform(
      matrix(1:10),
      matrix(10:1),
      NULL,
      TRUE,
      scale = NULL
    )
  )

  nu <- data.frame(x = 1:10)
  de <- data.frame(x = 10:1)
  ce <- data.frame(x = 1:5)

  d <- check.dataform(nu, de, scale = NULL, centers = ce, nullcenters = FALSE)

  expect_equal(d$nu, model.matrix(~., nu)[,-1, drop = FALSE], ignore_attr = TRUE)
  expect_equal(d$de, model.matrix(~., de)[,-1, drop = FALSE], ignore_attr = TRUE)
  expect_equal(d$ce, model.matrix(~., ce)[,-1, drop = FALSE], ignore_attr = TRUE)

  d <- check.dataform(nu, de, scale = "numerator", centers = ce, nullcenters = FALSE) |>
    suppressWarnings()

  nu_scale <- scale(nu)
  de_scale <- scale(de,
                    center = attr(nu_scale, "scaled:center"),
                    scale = attr(nu_scale, "scaled:scale"))
  ce_scale <- scale(ce,
                    center = attr(nu_scale, "scaled:center"),
                    scale = attr(nu_scale, "scaled:scale"))

  expect_equal(
    d$nu,
    model.matrix(~., as.data.frame(nu_scale))[,-1, drop = FALSE],
    ignore_attr = TRUE
  )
  expect_equal(
    d$de,
    model.matrix(~., as.data.frame(de_scale))[,-1, drop = FALSE],
    ignore_attr = TRUE
  )
  expect_equal(
    d$ce,
    model.matrix(~., as.data.frame(ce_scale))[,-1, drop = FALSE],
    ignore_attr = TRUE
  )

  d <- check.dataform(nu, de, scale = "denominator", centers = ce, nullcenters = FALSE) |>
    suppressWarnings()

  de_scale <- scale(de)
  nu_scale <- scale(nu,
                    center = attr(de_scale, "scaled:center"),
                    scale = attr(de_scale, "scaled:scale"))
  ce_scale <- scale(ce,
                    center = attr(de_scale, "scaled:center"),
                    scale = attr(de_scale, "scaled:scale"))

  expect_equal(
    d$nu,
    model.matrix(~., as.data.frame(nu_scale))[,-1, drop = FALSE],
    ignore_attr = TRUE
  )
  expect_equal(
    d$de,
    model.matrix(~., as.data.frame(de_scale))[,-1, drop = FALSE],
    ignore_attr = TRUE
  )
  expect_equal(
    d$ce,
    model.matrix(~., as.data.frame(ce_scale))[,-1, drop = FALSE],
    ignore_attr = TRUE
  )

  expect_warning(
    check.dataform(
      nu = data.frame(x = 1:10, y = 0),
      de = data.frame(x = 10:1, y = 1:10),
      scale = "numerator",
      centers = NULL,
      nullcenters = TRUE
    )
  )
  expect_warning(
    check.dataform(
      nu = data.frame(x = 1:10, y = 1:10),
      de = data.frame(x = 10:1, y = 1:10),
      scale = "numerator",
      centers = data.frame(x = 1:5, y = 1:5),
      nullcenters = FALSE
    )
  )
  expect_no_warning(
    check.dataform(
      nu = data.frame(x = 1:10, y = 0),
      de = data.frame(x = 10:1, y = 1:10),
      scale = NULL,
      centers = data.frame(x = 1:5, y = 1:5),
      nullcenters = FALSE
    )
  )
  expect_error(
    check.dataform(
      nu = data.frame(x = 1:10, y = 0),
      de = data.frame(x = 10:1, y = 1:10),
      scale = "both",
      centers = data.frame(x = 1:5, y = 1:5),
      nullcenters = FALSE
    )
  )

  expect_equal(
    check.dataform(
      nu = data.frame(x = 1:10, y = 1:10),
      de = data.frame(x = 10:1, y = 1:10),
      centers = NULL,
      nullcenters = TRUE,
      newdata = data.frame(x = 1:10, y = 1:10),
      scale = "numerator"
    ),
    scale(
      data.frame(x = 1:10, y = 1:10),
      center = c(5.5, 5.5),
      scale = c(sd(1:10), sd(1:10))
    ),
    ignore_attr = TRUE
  )
  expect_equal(
    check.dataform(
      nu = data.frame(x = 1:10, y = 1:10),
      de = data.frame(x = 10:1, y = 1:10),
      centers = NULL,
      nullcenters = TRUE,
      newdata = data.frame(x = 1:10, y = 1:10),
      scale = NULL
    ),
    data.frame(x = 1:10, y = 1:10) |> as.matrix(),
    ignore_attr = TRUE
  )

})

test_that("check.variables works", {
  expect_silent(
    check.variables(numerator_small, denominator_small)
  )
  expect_error(
    check.variables(numerator_small[,-1], denominator_small[,-2])
  )
  expect_error(
    check.variables(numerator_small[,-1], denominator_small[,-1], numerator_small[,-2])
  )
  expect_error(
    check.variables(numerator_small, cbind(as.factor(numerator_small$x1), numerator_small[,-1]))
  )
})

test_that("check.sigma works", {
  D1 <- distance(as.matrix(1:20), as.matrix(20:1))
  D2 <- distance(as.matrix(1:5), as.matrix(5:1))
  expect_type(
    check.sigma(10, NULL, NULL, D1),
    "double"
  )
  expect_warning(
    check.sigma(10, NULL, NULL, D2)
  )
  expect_warning(
    check.sigma(10, c(0.1, 0.1, 0.9), NULL, D1)
  )
  expect_error(
    check.sigma(10, NULL, c(-1, 1, 2), D2)
  )
  expect_error(
    check.sigma(10, NULL, c("a", "b"), D1)
  )
  expect_error(
    check.sigma(10, c(0,0.5,1), NULL, D2)
  )
  expect_error(
    check.sigma(10, c("a", "b"), NULL, D2)
  )
  expect_length(
    check.sigma(5, c(0.1,0.2), NULL, D1),
    2
  )
  expect_equal(
    check.sigma(5, c(0.5), NULL, D2),
    sqrt(median(D2[D2>0])/2)
  )
  expect_equal(
    check.sigma(1, NULL, NULL, D2),
    sqrt(median(D2[D2>0])/2)
  )
  expect_error(
    check.sigma("a", NULL, NULL, D2)
  )
  expect_error(
    check.sigma(0, NULL, NULL, D2)
  )
})

test_that("check.sigma_quantile.lhss works", {
  expect_no_error(
    check.sigma_quantile.lhss(10, NULL, NULL)
  )
  expect_true(
    all(check.sigma_quantile.lhss(10, NULL, NULL) > 0)
  )
  expect_true(
    all(check.sigma_quantile.lhss(10, NULL, NULL) < 1)
  )
  expect_null( # because sigma dominates everything, and this function only outputs probs for quantiles
    check.sigma_quantile.lhss(10, c(1,2,3), c(0.1, 0.2))
  )
  expect_type(
    check.sigma_quantile.lhss(10, NULL, NULL),
    "double"
  )
  expect_type(
    check.sigma_quantile.lhss(10, NULL, c(0.1, 0.2)),
    "double"
  )
  expect_equal(
    check.sigma_quantile.lhss(10, NULL, c(0.1, 0.2)),
    c(0.1, 0.2)
  )
  expect_error(
    check.sigma_quantile.lhss(10, NULL, 0)
  )
  expect_error(
    check.sigma_quantile.lhss(10, NULL, 1)
  )
  expect_equal(
    check.sigma_quantile.lhss(1, NULL, NULL),
    0.5
  )
  expect_error(
    check.sigma_quantile.lhss(10, matrix(1:4, 2), NULL)
  )
  expect_error(
    check.sigma_quantile.lhss(10, c(1,2,"a"), NULL)
  )
  expect_error(
    check.sigma_quantile.lhss(10, NULL, matrix(1:4, 2))
  )
  expect_error(
    check.sigma_quantile.lhss(10, NULL, c(1,2,"a"))
  )
  expect_error(
    check.sigma_quantile.lhss("a", NULL, NULL)
  )
  expect_error(
    check.sigma_quantile.lhss(c(1,2), NULL, NULL)
  )
  expect_error(
    check.sigma_quantile.lhss(0, NULL, NULL)
  )
})

test_that("check.lambda works", {
  expect_equal(
    check.lambda(10, NULL),
    10^seq(3, -3, length.out = 10)
  )
  expect_equal(
    check.lambda(10, c(1,2,3)),
    c(1,2,3)
  )
  expect_error(
    check.lambda(10, matrix(1:4, 2))
  )
  expect_error(
    check.lambda(10, c(1,2,"a"))
  )
  expect_error(
    check.lambda(c(1,2), NULL)
  )
})

test_that("check.centers works", {
  dat <- check.dataform(
    numerator_small,
    denominator_small,
    numerator_small,
    nullcenters = FALSE,
    scale = NULL
  )
  expect_equal(
    check.centers(dat$nu, dat$ce, 200),
    check.datatype(dat$nu),
    ignore_attr = TRUE
  )
  expect_error(
    check.centers(dat$nu, NULL, "a")
  )
  expect_error(
    check.centers(dat$nu, NULL, c(1,2))
  )
  expect_error(
    check.centers(dat$nu, NULL, -10)
  )
  expect_true(
    sum(
      duplicated(
        rbind(
          check.datatype(dat$nu),
          check.centers(dat$nu, NULL, 10)
        )
      )[51:60]
    ) == 10
  )
})

test_that("check.intercept works", {
  expect_type(
    check.intercept(TRUE),
    "logical"
  )
  expect_error(
    check.intercept("TRUE")
  )
  expect_error(
    check.intercept(NA)
  )
})

test_that("check.symmetric works", {
  dat <- check.dataform(
    numerator_small,
    denominator_small,
    numerator_small,
    nullcenters = FALSE,
    scale = NULL
  )
  expect_true(
    check.symmetric(
      dat$nu,
      dat$ce
    )
  )
  expect_false(
    check.symmetric(
      dat$nu,
      dat$de
    )
  )
})

test_that("check.parallel works", {
  p <- c(TRUE, FALSE, NA)
  nthreads <- c(1,2,10)
  iterator <- list(1, c(1,2,3))
  expect_true(
    check.parallel(p[1], nthreads[2], iterator[2][[1]])
  )
  expect_warning(
    expect_false(
      check.parallel(p[1], nthreads[1], iterator[2][[1]])
    )
  )
  expect_warning(
    check.parallel(p[1], nthreads[3], iterator[1][[1]])
  )
  expect_false(
    check.parallel(p[1], nthreads[3], iterator[1][[1]])
  ) |> suppressWarnings()
  expect_error(
    check.parallel(p[3], nthreads[3], iterator[2][[1]])
  )
})

test_that("check.threads works", {
  expect_warning(
    check.threads(FALSE, 10)
  )
  expect_equal(
    check.threads(FALSE, NULL),
    0
  )
  expect_equal(
    check.threads(TRUE, NULL),
    0
  )
  expect_equal(
    check.threads(TRUE, 10),
    10
  )
  expect_error(
    check.threads(TRUE, "a")
  )
  expect_equal(
    check.threads(TRUE, -1),
    1
  ) |> suppressWarnings()
  expect_error(
    check.threads(TRUE, c(1,2))
  )
  expect_warning(
    check.threads(TRUE, -1)
  )
})

test_that("check.epsilon works", {
  expect_equal(
    check.epsilon(0.1),
    0.1
  )
  expect_equal(
    check.epsilon(c(0.1, 0.2, 0.3)),
    c(0.1, 0.2, 0.3)
  )
  expect_error(
    check.epsilon(c(0.1, 0.2, -0.3))
  )
  expect_error(
    check.epsilon(c("a", "B"))
  )
  expect_error(
    check.epsilon(matrix(1:10))
  )
  expect_equal(
    check.epsilon(NULL),
    10^{1:-5}
  )
})

test_that("check.maxit works", {
  expect_equal(check.maxit(1000), 1000)
  expect_error(check.maxit("a"))
  expect_error(check.maxit(-1))
  expect_error(check.maxit(Inf))
  expect_error(check.maxit(c(1,2)))
})

test_that("check.nfold works", {
  expect_equal(check.nfold(FALSE, 5, 100), rep(0, 100))
  expect_equal(
    {set.seed(123); check.nfold(TRUE, 5, 100)},
    {set.seed(123); sample(rep_len(0:4, 100))}
  )
  expect_error(check.nfold(TRUE, "a", 100))
  expect_error(check.nfold(TRUE, c(1,2), 100))
  expect_error(check.nfold(TRUE, 1, 100))
  expect_error(check.nfold(TRUE, 101, 100))
})

test_that("check.sigma.predict works", {
  dr <- kliep(numerator_small, denominator_small, nsigma = 5)
  dr_nocv <- kliep(numerator_small, denominator_small, cv = FALSE)
  dr_onesigma <- kliep(numerator_small, denominator_small, cv = FALSE, sigma = 1)
  expect_equal(
    check.sigma.predict(dr, 10),
    10
  )
  expect_equal(
    check.sigma.predict(dr, c(1,2,3)),
    c(1,2,3)
  )
  expect_equal(
    check.sigma.predict(dr, "all"),
    dr$sigma
  )
  expect_equal(
    check.sigma.predict(dr, "sigmaopt"),
    dr$sigma_opt
  )
  expect_error(
    check.sigma.predict(dr, c(1,2,"a"))
  )
  expect_error(
    check.sigma.predict(dr, "b")
  )
  expect_warning(
    check.sigma.predict(dr_nocv, "sigmaopt")
  )
  expect_equal(
    check.sigma.predict(dr_nocv, "sigmaopt"),
    dr_nocv$sigma
  ) |> suppressWarnings()
  expect_equal(
    check.sigma.predict(dr_onesigma, "sigmaopt"),
    dr_onesigma$sigma
  )
  expect_error(
    check.sigma.predict(dr, matrix(1:3))
  )
})

test_that("check.lambdasigma.predict works", {
  dr <- lhss(numerator_small, denominator_small, nsigma = 5, nlambda = 5)

  expect_equal(
    check.lambdasigma.predict(dr, "sigmaopt", c(1, 2), lambdaind = match(c(1,2), dr$lambda))[,4],
    c(dr$sigma[which.min(dr$cv_score[, which(dr$lambda == 1)]),
               which(dr$lambda == 1)],
      NA)
  )

  expect_equal(
    check.lambdasigma.predict(dr, "sigmaopt", dr$lambda_opt, match(dr$lambda_opt, dr$lambda))[,3:4],
    c(dr$lambda_opt, dr$sigma_opt),
    ignore_attr = TRUE
  )

  expect_equal(
    check.lambdasigma.predict(dr, "all", dr$lambda_opt, match(dr$lambda_opt, dr$lambda))[,3:4],
    matrix(
      c(rep(dr$lambda_opt, nrow(dr$sigma)),
        dr$sigma[, which(dr$lambda == dr$lambda_opt)]
       ),
      ncol = 2
    ),
    ignore_attr = TRUE
  )
  expect_equal(
    check.lambdasigma.predict(dr, "all", c(1,2), match(c(1,2), dr$lambda))[,3:4],
    matrix(
      c(rep(c(1,2), each = nrow(dr$sigma)),
        dr$sigma[, match(c(1,2), dr$lambda)]
       ),
      ncol = 2
    ),
    ignore_attr = TRUE
  )
  expect_equal(
    check.lambdasigma.predict(
      dr, c(1,2), c(1,2), match(c(1,2), dr$lambda)
    ),
    matrix(c(3, 3, rep(NA, 6), 1,1,2,2,1,2,1,2), ncol = 4),
    ignore_attr = TRUE
  )
  expect_error(
    check.lambdasigma.predict(dr, "b", c(1,2), lambdaind = match(c(1,2), dr$lambda))
  )
  expect_error(
    check.lambdasigma.predict(dr, matrix(1:3), c(1,2), lambdaind = match(c(1,2), dr$lambda))
  )
})

test_that("check.lambda.predict works", {
  dr <- ulsif(numerator_small, denominator_small, nlambda = 5)
  expect_equal(
    check.lambda.predict(dr, "lambdaopt"),
    dr$lambda_opt
  )
  expect_equal(
    check.lambda.predict(dr, "all"),
    dr$lambda
  )
  expect_equal(
    check.lambda.predict(dr, c(1,2,3)),
    c(1,2,3)
  )
  expect_error(
    check.lambda.predict(dr, c(1,2,"a"))
  )
  expect_error(
    check.lambda.predict(dr, matrix(1:10))
  )
})

test_that("check.subspace.spectral.predict works", {
  dr <- spectral(numerator_small, denominator_small, m = 1:10)
  expect_equal(
    check.subspace.spectral.predict(dr, "opt"),
    dr$m_opt
  )
  expect_equal(
    check.subspace.spectral.predict(dr, "all"),
    dr$m
  )
  expect_equal(
    check.subspace.spectral.predict(dr, c(1,2,3)),
    c(1,2,3)
  )
  expect_error(
    check.subspace.spectral.predict(dr, c(1,2,101))
  )
  expect_error(
    check.subspace.spectral.predict(dr, c(1,2,"a"))
  )
  expect_error(
    check.subspace.spectral.predict(dr, (matrix(10, )))
  )
})

test_that("check.subspace works", {
  expect_equal(
    check.subspace(10, 100),
    10
  )
  expect_equal(
    check.subspace(NULL, 100),
    10
  )
  expect_error(
    check.subspace("a", 100)
  )
  expect_error(
    check.subspace(1.1, 100)
  )
  expect_error(
    check.subspace(11, 10)
  )
})

test_that("check.subspace.spectral works", {
  expect_equal(
    check.subspace.spectral(1:99, 1:100),
    1:99
  )
  expect_equal(
    check.subspace.spectral(NULL, rep(1:10, each = 10)),
    unique(floor(seq(1, 90, length.out = 50)))
  )
  expect_error(
    check.subspace.spectral(1:100, rep(1:10, each = 10))
  )
  expect_error(
    check.subspace.spectral(c(-1, 10), 1:100)
  )
  expect_error(
    check.subspace.spectral(c("a", "b"), 1:100)
  )
})

test_that("check.newdata works", {
  dr <- kliep(numerator_small, denominator_small)

  expect_equal(
    check.newdata(dr, numerator_small),
    dr$model_matrices$nu,
    ignore_attr = TRUE
  )
  expect_error(
    check.newdata(dr, denominator_small[,c(1,3,2)])
  )
})

test_that("check.var.names works", {
  expect_silent(
    check.var.names(c("x1", "x2"), numerator_small)
  )
  expect_error(
    check.var.names(c("x1", "X2"), denominator_small)
  )
})

test_that("check.object.type works", {
  expect_silent(
    check.object.type(
      ulsif(numerator_small, denominator_small, nsigma = 5, nlambda = 5)
    )
  )
  expect_error(
    check.object.type(
      data.frame(x = 1, y = 2)
    )
  )
})

test_that("check.logscale works", {
  ext <- data.frame(dr = c(-0.01, 1, 2, 1))

  expect_warning(
    check.logscale(ext, TRUE, tol = 1e-6)
  )
  expect_type(
    check.logscale(ext, FALSE, tol = 1e-6),
    "list"
  )
})

