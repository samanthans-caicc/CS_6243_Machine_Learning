# -----------------------------------------------------------------------------
# PROBLEM 4 - Univariate polynomial regression on each Auto-MPG feature
# -----------------------------------------------------------------------------
# Requires, in order: setup.R, q1.R (standardize), q2.R (poly_features,
# standardize_poly, solve_LS, mse, fit_poly), q3.R (TRI_raw, TRO, TEI_raw,
# TEO, feature_names, cor_matrix).

library(ggplot2)
library(tidyr)

degrees <- 0:5

# predict_poly: evaluate a fitted model at NEW raw x values, using only the
# TRAINING statistics stored in the fit (x mean/sd, then column mean/sd).
predict_poly <- function(fit, x_new_raw, d) {
  z <- (x_new_raw - fit$x_stats$mu) / fit$x_stats$sigma     # training x stats
  X <- poly_features(z, d)
  ps <- fit$poly_stats
  cols <- which(ps$standardize_column)
  if (length(cols) > 0) {                                     # training column stats
    X[, cols] <- sweep(sweep(X[, cols, drop = FALSE], 2, ps$mu[cols], "-"),
                       2, ps$sigma[cols], "/")
  }
  as.numeric(X %*% fit$w)
}

# 1. Double loop over features and degrees.
#    fit_poly is given the RAW feature row; it standardizes x with training
#    stats itself, so the policy is identical to Problem 2.
train_mse_mat <- matrix(NA_real_, nrow = 7, ncol = length(degrees),
                        dimnames = list(feature_names, paste0("M=", degrees)))
test_mse_mat <- train_mse_mat
fits <- list()   # keep every fit so the figures can reuse them

for (j in seq_along(feature_names)) {
  for (M in degrees) {
    fit <- fit_poly(TRI_raw[j, ], TRO[, 1], TEI_raw[j, ], TEO[, 1], M)
    train_mse_mat[j, M + 1] <- fit$train_mse
    test_mse_mat[j, M + 1]  <- fit$test_mse
    fits[[paste(feature_names[j], M)]] <- fit
  }
}

# Sanity check: M = 0 is the bias-only model (predicts mean of TRO),
# so every feature must give the same MSE.
stopifnot(diff(range(test_mse_mat[, "M=0"])) < 1e-10)

cat("Training MSE (mpg^2):\n"); print(round(train_mse_mat, 3))
cat("Test MSE (mpg^2):\n");     print(round(test_mse_mat, 3))

best <- which(test_mse_mat == min(test_mse_mat), arr.ind = TRUE)[1, ]
cat(sprintf("Best univariate model: %s, M = %d, test MSE = %.3f\n",
            feature_names[best[1]], degrees[best[2]], min(test_mse_mat)))

# Compare with Problem 3 correlations (Task 4.4 discussion)
comparison <- data.frame(abs_corr_with_mpg = round(cor_matrix["mpg", feature_names], 3),
                         best_test_mse     = round(apply(test_mse_mat, 1, min), 3))
print(comparison[order(-comparison$abs_corr_with_mpg), ])

# Helper for Figs. 8 and 9: fitted curve on a dense RAW grid + training scatter.
# The scatter uses RAW feature values so it shares the x-axis with the curve.
plot_fit_curve <- function(feature, M, grid_range, x_label, fig_title, file) {
  j <- match(feature, feature_names)
  fit <- fits[[paste(feature, M)]]
  grid <- seq(grid_range[1], grid_range[2], length.out = 500)
  curve_df <- data.frame(x = grid, mpg = predict_poly(fit, grid, M))
  scatter_df <- data.frame(x = TRI_raw[j, ], mpg = TRO[, 1])

  p <- ggplot() +
    geom_point(data = scatter_df, aes(x = x, y = mpg, colour = "Training data"),
               alpha = 0.6) +
    geom_line(data = curve_df, aes(x = x, y = mpg,
                                   colour = sprintf("Degree-%d fit", M)),
              linewidth = 1) +
    scale_colour_manual(values = setNames(c("blue", "red"),
                                          c("Training data", sprintf("Degree-%d fit", M)))) +
    labs(x = x_label, y = "Fuel economy (mpg)", colour = NULL, title = fig_title) +
    theme(legend.position = "top")
  ggsave(file, p, width = 6.5, height = 4.5, dpi = 150)
  p
}

# Task 4.1: Fig. 8 - weight, M = 3, grid over [1613, 5140]
fig8 <- plot_fit_curve("weight", 3, c(1613, 5140), "Weight (lb)",
                       "Fig. 8: mpg vs weight, cubic fit (M = 3)",
                       "fig8_weight_M3.png")

# Task 4.2: Fig. 9 - horsepower, M = 3, grid over [46, 230]
fig9 <- plot_fit_curve("horsepower", 3, c(46, 230), "Horsepower (hp)",
                       "Fig. 9: mpg vs horsepower, cubic fit (M = 3)",
                       "fig9_horsepower_M3.png")

# Task 4.3: Fig. 10 - test MSE vs M, one coloured line per feature
mse_long <- data.frame(feature = factor(rep(feature_names, times = length(degrees)),
                                        levels = feature_names),
                       M = rep(degrees, each = 7),
                       test_mse = as.vector(test_mse_mat))
fig10 <- ggplot(mse_long, aes(x = M, y = test_mse, colour = feature)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = degrees) +
  labs(x = "Polynomial degree M", y = "Test MSE (mpg^2)", colour = "Feature",
       title = "Fig. 10: Test MSE vs degree for each single feature")
ggsave("fig10_test_mse_vs_M.png", fig10, width = 7, height = 4.5, dpi = 150)

# Task 4.4: Discussion goes in the report (use `comparison` and Fig. 10).