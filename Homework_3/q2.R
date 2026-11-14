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
# as.numeric() on both: X %*% w is an n x 1 matrix, and y from numpy is a
# 1-D array, which R refuses to subtract from a matrix.
mse <- function(yhat, y) {
  mean((as.numeric(yhat) - as.numeric(y))^2)
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
library(ggplot2)

# Task 2.1: Fig. 1 - standardized x_train vs y_train (train_size = 5000)
split_5000 <- split_data(5000)
z_5000 <- standardize(split_5000$x_train, split_5000$x_test)
fig1 <- ggplot(data.frame(x = z_5000$train, y = split_5000$y_train), aes(x = x, y = y)) +
  geom_point(alpha = 0.3, size = 0.8) +
  labs(x = "MedInc (standardized, unitless)",
       y = "Median house value ($100k)",
       title = "Fig. 1: Training data, standardized MedInc vs house value (train_size = 5000)")
ggsave("fig1_standardized_scatter.png", fig1, width = 6.5, height = 4.5, dpi = 150)

# Task 2.2: train/test MSE vs degree d (train_size = 5000)
d_values <- c(1, 2, 3, 4, 5, 10, 15, 20)
mse_vs_d <- data.frame(d = d_values, t(sapply(d_values, function(d) {
  fit <- fit_poly(split_5000$x_train, split_5000$y_train,
                  split_5000$x_test, split_5000$y_test, d)
  c(train_mse = fit$train_mse, test_mse = fit$test_mse)
})))
cat("Problem 2 - MSE vs degree d (train_size = 5000), units ($100k)^2:\n")
print(transform(mse_vs_d, train_mse = round(train_mse, 4), test_mse = round(test_mse, 4)))

# Task 2.3: Fig. 2 - train (blue) and test (red) MSE vs d
mse_d_long <- data.frame(d = rep(d_values, 2),
                         set = rep(c("Train", "Test"), each = length(d_values)),
                         mse = c(mse_vs_d$train_mse, mse_vs_d$test_mse))
fig2 <- ggplot(mse_d_long, aes(x = d, y = mse, colour = set)) +
  geom_line() +
  geom_point() +
  scale_colour_manual(values = c(Train = "blue", Test = "red")) +
  scale_x_continuous(breaks = d_values) +
  labs(x = "Polynomial degree d", y = "MSE (($100k)^2)", colour = NULL,
       title = "Fig. 2: Train and test MSE vs degree (train_size = 5000)")
ggsave("fig2_mse_vs_degree.png", fig2, width = 6.5, height = 4.5, dpi = 150)

# Task 2.4: d = 9, nested training sizes. fit_poly recomputes all training
#           statistics (mu_x, sigma_x, mu_k, sigma_k) for each split.
n_values <- c(1000, 3000, 5000, 7000, 10000)
mse_vs_n <- data.frame(train_size = n_values, t(sapply(n_values, function(n) {
  s <- split_data(n)
  fit <- fit_poly(s$x_train, s$y_train, s$x_test, s$y_test, 9)
  c(train_mse = fit$train_mse, test_mse = fit$test_mse)
})))
cat("\nProblem 2 - MSE vs training size (d = 9), units ($100k)^2:\n")
print(transform(mse_vs_n, train_mse = round(train_mse, 4), test_mse = round(test_mse, 4)))

# Task 2.5: Fig. 3 - test MSE vs training size (d = 9)
fig3 <- ggplot(mse_vs_n, aes(x = train_size, y = test_mse)) +
  geom_line(colour = "red") +
  geom_point(colour = "red") +
  scale_x_continuous(breaks = n_values) +
  labs(x = "Training set size (samples)", y = "Test MSE (($100k)^2)",
       title = "Fig. 3: Test MSE vs training size (d = 9)")
ggsave("fig3_test_mse_vs_n.png", fig3, width = 6.5, height = 4.5, dpi = 150)

# Task 2.6: Fig. 4 - d = 9 fit over the training scatter, sorted by x_train.
#           Fitted values are X_poly w on the training set itself.
fit_d9 <- fit_poly(split_5000$x_train, split_5000$y_train,
                   split_5000$x_test, split_5000$y_test, 9)
order_x <- order(z_5000$train)
curve_df <- data.frame(x = z_5000$train[order_x],
                       y = as.numeric(fit_d9$poly_stats$X_train %*% fit_d9$w)[order_x])
fig4 <- ggplot() +
  geom_point(data = data.frame(x = z_5000$train, y = split_5000$y_train),
             aes(x = x, y = y, colour = "Training data"), alpha = 0.3, size = 0.8) +
  geom_line(data = curve_df, aes(x = x, y = y, colour = "Degree-9 fit"), linewidth = 1) +
  scale_colour_manual(values = c("Training data" = "blue", "Degree-9 fit" = "red")) +
  labs(x = "MedInc (standardized, unitless)", y = "Median house value ($100k)",
       colour = NULL, title = "Fig. 4: Degree-9 polynomial fit (train_size = 5000)") +
  theme(legend.position = "top")
ggsave("fig4_degree9_fit.png", fig4, width = 6.5, height = 4.5, dpi = 150)

# Task 2.7: Discussion of Figs. 1-4 goes in the report.
