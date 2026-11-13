# -----------------------------------------------------------------------------
# BONUS 2 - KNN regression (optional)
# -----------------------------------------------------------------------------
# Requires, in order: setup.R, q1.R (split_data, standardize), q2.R (mse).
# Optional: bq1.R first, so the comparison with radius regression prints.
# Implemented in base R (the HW says NumPy; state this in the report).

library(ggplot2)

# x_train, x_test are plain numeric VECTORS (one scalar feature), so the
# Euclidean distance is just |x_j - x_i|.
# Ties: order() is stable, so equally distant points are taken in training-
# index order (state this tie-breaking rule in the report).
knn_predict <- function(x_train, y_train, x_test, K) {
  K_actual <- min(K, length(x_train))      # K > n_train -> use all points
  predictions <- numeric(length(x_test))
  for (i in seq_along(x_test)) {
    distances <- abs(x_train - x_test[i])
    neighbors <- order(distances)[seq_len(K_actual)]
    predictions[i] <- mean(y_train[neighbors])
  }
  list(predictions = predictions)
}

# Split + standardize with TRAINING stats, then run KNN
run_knn <- function(train_size, K) {
  s <- split_data(train_size)
  z <- standardize(s$x_train, s$x_test)
  p <- knn_predict(z$train, s$y_train, z$test, K)$predictions
  mse(p, s$y_test)
}

# Task B2.1: train_size = 5000, sweep K
K_values <- c(1, 3, 5, 10, 20, 50, 100, 500, 1000)
knn_vs_K <- data.frame(K = K_values,
                       test_mse = sapply(K_values, function(K) run_knn(5000, K)))
cat("Bonus 2 - test MSE vs K (train_size = 5000):\n")
print(transform(knn_vs_K, test_mse = round(test_mse, 4)))

# Task B2.2: Fig. 14 - test MSE vs K (log x-axis: K spans 1 to 1000)
fig14 <- ggplot(knn_vs_K, aes(x = K, y = test_mse)) +
  geom_line(colour = "blue") +
  geom_point(colour = "blue") +
  scale_x_log10(breaks = K_values) +
  labs(x = "Number of neighbours K (log scale)",
       y = "Test MSE (($100k)^2)",
       title = "Fig. 14: KNN regression, test MSE vs K (train_size = 5000)")
ggsave("fig14_knn_mse_vs_K.png", fig14, width = 7, height = 4.5, dpi = 150)

# Task B2.3: K = 100, nested training prefixes, fixed test set
train_sizes <- seq(1000, 10000, by = 1000)
knn_vs_n <- data.frame(train_size = train_sizes,
                       test_mse = sapply(train_sizes, function(n) run_knn(n, 100)))
cat("\nBonus 2 - test MSE vs training size (K = 100):\n")
print(transform(knn_vs_n, test_mse = round(test_mse, 4)))

# Task B2.4: Fig. 15 - test MSE vs training size
fig15 <- ggplot(knn_vs_n, aes(x = train_size, y = test_mse)) +
  geom_line(colour = "blue") +
  geom_point(colour = "blue") +
  scale_x_continuous(breaks = train_sizes) +
  labs(x = "Training set size (samples)",
       y = "Test MSE (($100k)^2)",
       title = "Fig. 15: KNN regression, test MSE vs training size (K = 100)")
ggsave("fig15_knn_mse_vs_n.png", fig15, width = 7, height = 4.5, dpi = 150)

# Task B2.5: numbers for the comparison (discussion goes in the report)
best_knn <- knn_vs_K[which.min(knn_vs_K$test_mse), ]
cat(sprintf("\nBest KNN (train_size = 5000): K = %d, test MSE = %.4f\n",
            best_knn$K, best_knn$test_mse))
if (exists("radius_vs_C")) {
  best_rad <- radius_vs_C[which.min(radius_vs_C$test_mse), ]
  cat(sprintf("Best radius (train_size = 5000): C = %g, test MSE = %.4f\n",
              best_rad$C, best_rad$test_mse))
}