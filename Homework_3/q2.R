# -----------------------------------------------------------------------------
# PROBLEM 2 - Polynomial regression, one feature (Dataset A)
# -----------------------------------------------------------------------------
# poly_features(x, d)
#   - columns x^0, x^1, ..., x^d
#   - note in the report which orientation you used (HW writes (d+1) x n)
poly_features <- function(x, d) {
  sapply(0:d, function(power) x^power)
  # Return the matrix with columns x^0, x^1, ..., x^d
}
 
# Column (row, in HW notation) standardizer
#   - skip the bias column
#   - skip any column with sd = 0
#   - learn stats on train, apply to test, return the stats
standardize_poly <- function(X_train, X_test) {
  X_train <- as.matrix(X_train)
  X_test <- as.matrix(X_test)
  n <- nrow(X_train)
  mu <- colMeans(X_train)
  sigma <- apply(X_train, 2, function(column) sqrt(mean((column - mu[col(column)])^2)))
  standardize_column <- sigma > 0

  X_train <- X_train
  X_test <- X_test
  X_train[, standardize_column] <- (X_train[, standardize_column] - mu[standardize_column]) / sigma[standardize_column]
  X_test[, standardize_column] <- (X_test[, standardize_column] - mu[standardize_column]) / sigma[standardize_column]

  list(
    X_train = X_train,
    X_test = X_test,
    mu = mu,
    sigma = sigma,
    standardize_column = standardize_column
  )
}
 
# solve_LS(X, y) via SVD
#   - take the SVD of X
#   - drop singular values below a tolerance (report the tolerance you use)
#   - w = V * D^-1 * U^T * y
solve_LS <- function(X, y) {
  svd_result <- svd(X)
  tolerance <- max(dim(X)) * max(svd_result$d) * .Machine$double.eps
  keep <- svd_result$d > tolerance
  svd_result$u[, keep, drop = FALSE] %*%
    (svd_result$v[, keep, drop = FALSE] %*%
     (t(svd_result$u[, keep, drop = FALSE]) %*% y / svd_result$d[keep]))
}
 
# mse(yhat, y) - mean squared error in ORIGINAL y units
mse <- function(yhat, y) {
  mean((yhat - y)^2)

}
 
# Core pipeline (reused in Problems 2, 4, 5):
#   given a train/test split and degree d:
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
 
# Task 2.1: Fig. 1 - standardized x_train vs y_train (train_size = 5000).
# Task 2.2: train_size = 5000, loop d in {1,2,3,4,5,10,15,20}; table of
#           train/test MSE.
# Task 2.3: Fig. 2 - train MSE (blue) and test MSE (red) vs d.
# Task 2.4: d = 9, loop train_size in {1000,3000,5000,7000,10000};
#           RE-SPLIT each time so all stats are recomputed. Table.
# Task 2.5: Fig. 3 - test MSE vs train_size (d = 9).
# Task 2.6: Fig. 4 - d = 9, train_size = 5000: scatter (blue) + fitted
#           curve (red). Sort by x_train before drawing the line.
# Task 2.7: Discussion - data shape, MSE vs d, MSE vs train size,
#           under/overfitting, quality of the d = 9 fit.
 
 