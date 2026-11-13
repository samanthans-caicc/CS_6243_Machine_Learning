# -----------------------------------------------------------------------------
# BONUS 1 - Radius-based regression (optional)
# -----------------------------------------------------------------------------
# Requires, in order: setup.R, q1.R (split_data, standardize), q2.R (mse).
# Implemented in base R (the HW says NumPy; state in the report that the
# same algorithm was written from scratch in R instead).

library(ggplot2)

# x_train, x_test are plain numeric VECTORS (one scalar feature, MedInc),
# so the distance is just |x_j - x_i| - no matrix rows needed.
radius_predict <- function(x_train, y_train, x_test, C) {
  predictions <- numeric(length(x_test))
  used_fallback <- logical(length(x_test))
  y_mean <- mean(y_train)

  for (i in seq_along(x_test)) {
    neighbors <- abs(x_train - x_test[i]) <= C
    if (any(neighbors)) {
      predictions[i] <- mean(y_train[neighbors])
    } else {
      predictions[i] <- y_mean          # explicit fallback
      used_fallback[i] <- TRUE
    }
  }
  list(predictions = predictions,
       fallback_count = sum(used_fallback),
       fallback_fraction = mean(used_fallback))
}

# Split + standardize with TRAINING stats, then run radius regression
run_radius <- function(train_size, C) {
  s <- split_data(train_size)
  z <- standardize(s$x_train, s$x_test)
  r <- radius_predict(z$train, s$y_train, z$test, C)
  c(test_mse = mse(r$predictions, s$y_test),
    fallback_fraction = r$fallback_fraction)
}

# Task B1.1: train_size = 5000, sweep C (standardized x-scale)
C_values <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.3, 0.5, 1, 2, 3, 4, 5)
radius_vs_C <- data.frame(C = C_values,
                          t(sapply(C_values, function(C) run_radius(5000, C))))
cat("Bonus 1 - test MSE vs C (train_size = 5000):\n")
print(transform(radius_vs_C, test_mse = round(test_mse, 4),
                             fallback_fraction = round(fallback_fraction, 4)))

# Task B1.2: Fig. 12 - test MSE vs C (log x-axis: C spans 0.001 to 5)
fig12 <- ggplot(radius_vs_C, aes(x = C, y = test_mse)) +
  geom_line(colour = "blue") +
  geom_point(colour = "blue") +
  scale_x_log10(breaks = C_values, labels = as.character(C_values)) +
  labs(x = "Radius C (standardized MedInc units, log scale)",
       y = "Test MSE (($100k)^2)",
       title = "Fig. 12: Radius regression, test MSE vs C (train_size = 5000)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("fig12_radius_mse_vs_C.png", fig12, width = 7, height = 4.5, dpi = 150)

# Task B1.3: C = 0.1, nested training prefixes, fixed test set.
#            standardize() is recomputed for every train_size inside run_radius.
train_sizes <- seq(1000, 10000, by = 1000)
radius_vs_n <- data.frame(train_size = train_sizes,
                          t(sapply(train_sizes, function(n) run_radius(n, 0.1))))
cat("\nBonus 1 - test MSE vs training size (C = 0.1):\n")
print(transform(radius_vs_n, test_mse = round(test_mse, 4),
                             fallback_fraction = round(fallback_fraction, 4)))

# Task B1.4: Fig. 13 - test MSE vs training size
fig13 <- ggplot(radius_vs_n, aes(x = train_size, y = test_mse)) +
  geom_line(colour = "blue") +
  geom_point(colour = "blue") +
  scale_x_continuous(breaks = train_sizes) +
  labs(x = "Training set size (samples)",
       y = "Test MSE (($100k)^2)",
       title = "Fig. 13: Radius regression, test MSE vs training size (C = 0.1)")
ggsave("fig13_radius_mse_vs_n.png", fig13, width = 7, height = 4.5, dpi = 150)

# Task B1.5: the fallback fractions are the `fallback_fraction` columns above;
#            the discussion goes in the report.