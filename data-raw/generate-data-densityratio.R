## code to prepare `numerator_data` and `denominator_data` datasets


set.seed(1)

library(tibble)
library(dplyr)
library(ggplot2)

n <- 1000
numerator_data <- tibble(
  x1 = sample(c("A", "B", "C"), n, replace = TRUE) |> factor(),
  x2 = rbinom(n, 1, prob = 1 / (1 + exp(-(as.numeric(x1) - mean(as.numeric(x1)))/2))) |>
    factor(labels = c("G1", "G2")),
  x3means = model.matrix( ~ interaction(x1, x2)) %*% c(2:7/10) |> c(),
  x3 = x3means + rnorm(n, 0, sqrt(1 - var(x3means))),
  x4 = 0.3 * x3 + rnorm(n, 0, sqrt(1 - 0.3^2)),
  x5group = rbinom(n, 1, 0.4),
  x5 = x5group * rnorm(n, 0, 1) + (1 - x5group) * rnorm(n, 3, 1)
) |> select(-c(x3means, x5group))

denominator_data <- tibble(
  x1 = sample(c("A", "B", "C"), n, replace = TRUE, prob = c(0.25, 0.25, 0.5)) |> factor(),
  x2 = rbinom(n, 1, prob = 0.5) |> factor(labels = c("G1", "G2")),
  x3means = model.matrix(~interaction(x1, x2)) %*% c(2:7/10) |> c(),
  x3 = -0.5 + x3means + rnorm(n, 0, sqrt(2 - var(x3means))),
  x4 = rnorm(n, 1, 1),
  x5 = rnorm(n, 3, 2)
) |> select(-x3means)

usethis::use_data(numerator_data, denominator_data, overwrite = TRUE)
