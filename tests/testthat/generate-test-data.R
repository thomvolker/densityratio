# create some test data
set.seed(45)
N <- 200
P_max <- 100
S_denominator <- rWishart(1, P_max, diag(P_max))[,,1]
test_df_numerator_max <- data.frame(matrix(rnorm(N*P_max), N))
test_df_denominator_max <- data.frame(matrix(rnorm(N*P_max), N) %*% chol(S_denominator))
test_df_numerator_10   <- test_df_numerator_max[, 1:10]
test_df_denominator_10 <- test_df_denominator_max[, 1:10]
test_df_numerator_1    <- test_df_numerator_max[, 1, drop = FALSE]
test_df_denominator_1  <- test_df_denominator_max[, 1, drop = FALSE]

# get some test values in here
# the density ratio at 0 should be greater than the density ratio at -10
# for manual checking, also include colmeans at numerator and denominator
test_df_newdata_max <- data.frame(rbind(
  t(rep(0, P_max)),
  t(rep(-10, P_max)),
  t(colMeans(test_df_numerator_max)),
  t(colMeans(test_df_denominator_max))
))
test_df_newdata_10 <- test_df_newdata_max[, 1:10]
test_df_newdata_1 <- test_df_newdata_max[, 1, drop = FALSE]
