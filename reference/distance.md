# Create a Gram matrix with squared Euclidean distances between observations in the input matrix `X` and the input matrix `Y`

Create a Gram matrix with squared Euclidean distances between
observations in the input matrix `X` and the input matrix `Y`

## Arguments

- X:

  A numeric input matrix

- Y:

  A numeric input matrix with the same variables as `X`

- intercept:

  Logical indicating whether an intercept should be added to the
  estimation procedure. In this case, the first column is an all-zero
  column (which will be transformed into an all-ones column in the
  kernel).
