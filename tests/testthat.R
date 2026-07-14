# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)

Sys.setenv(OMP_THREAD_LIMIT = "1", OPENBLAS_NUM_THREADS = "1")
if(requireNamespace("RcppArmadillo", quietly = TRUE)) {
  RcppArmadillo::armadillo_throttle_cores()
}
library(densityratio)
test_check("densityratio")
if(requireNamespace("RcppArmadillo", quietly = TRUE)) {
  RcppArmadillo::armadillo_reset_cores()
}
Sys.unsetenv(c("OMP_THREAD_LIMIT", "OPENBLAS_NUM_THREADS"))
