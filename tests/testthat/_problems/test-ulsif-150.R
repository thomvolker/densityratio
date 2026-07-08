# Extracted from test-ulsif.R:150

# test -------------------------------------------------------------------------
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
