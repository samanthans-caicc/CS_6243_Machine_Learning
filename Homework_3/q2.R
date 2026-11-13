# -----------------------------------------------------------------------------
# PROBLEM 2 - Polynomial regression, one feature (Dataset A)
# -----------------------------------------------------------------------------
# poly_features(x, d)
#   - columns x^0, x^1, ..., x^d
#   - note in the report which orientation you used (HW writes (d+1) x n)
# Orientation used: n x (d+1) (one row per sample) = transpose of the HW's
# (d+1) x n. outer() always returns a matrix, even when d = 0 or n = 1.
poly_features <- function(x, d) {
  outer(as.numeric(x), 0:d, "^")
}

# Column (row, in HW notation) standardizer
#   - skip the bias column
#   - skip any column with sd = 0
#   - learn stats on train, apply to test, return the stats
standardize_poly <- function(X_train, X_test) {
  X_train <- as.matrix(X_train)
  X_test <- as.matrix(X_test)
  mu <- colMeans(X_train)
  # Population sd (divide by n) of each column
  sigma <- apply(X_train, 2, function(column) sqrt(mean((column - mean(column))^2)))

  # Skip the bias column (column 1) explicitly, and any constant column (sd = 0)
  standardize_column <- sigma > 0
  standardize_column[1] <- FALSE
  cols <- which(standardize_column)

  # sweep() subtracts/divides each column by ITS OWN mu/sigma.
  # (Plain X - mu would recycle mu down the rows, which is wrong.)
  if (length(cols) > 0) {
    X_train[, cols] <- sweep(sweep(X_train[, cols, drop = FALSE], 2, mu[cols], "-"),
                             2, sigma[cols], "/")
    X_test[, cols]  <- sweep(sweep(X_test[, cols, drop = FALSE], 2, mu[cols], "-"),
                             2, sigma[cols], "/")
  }

  list(
    X_train = X_train,
    X_test = X_test,
    mu = mu,
    sigma = sigma,
    standardize_column = standardize_column
  )
}

# solve_LS(X, y) via SVD
# X = U D V^T  =>  w = V D^-1 U^T y  (pseudo-inverse solution)
# Tolerance: singular values <= max(n, p) * max(d) * eps are treated as zero
# (same default rule as numpy.linalg.pinv / MATLAB pinv).
solve_LS <- function(X, y) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  s <- svd(X)
  tolerance <- max(dim(X)) * max(s$d) * .Machine$double.eps
  keep <- s$d > tolerance
  U <- s$u[, keep, drop = FALSE]
  V <- s$v[, keep, drop = FALSE]
  V %*% ((t(U) %*% y) / s$d[keep])   # (p x 1) weight vector
}

# mse(yhat, y) - mean squared error in ORIGINAL y units
mse <- function(yhat, y) {
  mean((yhat - y)^2)
}

# Core pipeline (reused in Problems 2, 4, 5):
#   standardize x -> build poly features -> standardize columns
#   -> solve LS -> return train MSE, test MSE, w (and stats if needed)
fit_poly <- function(x_train, y_train, x_test, y_test, d) {
  x_stats <- standardize(x_train, x_test)
  X_train <- poly_features(x_stats$train, d)
  X_test <- poly_features(x_stats$test, d)
  poly_stats <- standardize_poly(X_train, X_test)
  w <- solve_LS(poly_stats$X_train, y_train)

  list(
    train_mse = mse(poly_stats$X_train %*% w, y_train),
    test_mse = mse(poly_stats$X_test %*% w, y_test),
    w = w,
    x_stats = x_stats,
    poly_stats = poly_stats
  )
}