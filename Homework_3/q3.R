# -----------------------------------------------------------------------------
# PROBLEM 3 - Auto-MPG standardization and correlation (Dataset B)
# -----------------------------------------------------------------------------
# Requires setup.R to have been sourced first (uses `auto` and `perm_auto`).

library(ggplot2)
library(tidyr)

feature_names <- c("cylinders", "displacement", "horsepower", "weight",
                   "acceleration", "model_year", "origin")

# 1. Split using the numpy seed-42 permutation from setup.R
train_rows <- perm_auto[1:300]
test_rows  <- perm_auto[301:392]

# 2. Arrange as TRI, TRO, TEI, TEO (features in ROWS, as in the HW)
TRI <- t(as.matrix(auto[train_rows, feature_names]))   # 7 x 300
TRO <- matrix(auto$mpg[train_rows], ncol = 1)          # 300 x 1
TEI <- t(as.matrix(auto[test_rows, feature_names]))    # 7 x 92
TEO <- matrix(auto$mpg[test_rows], ncol = 1)           # 92 x 1
colnames(TRI) <- NULL
colnames(TEI) <- NULL
stopifnot(all(dim(TRI) == c(7, 300)), all(dim(TEI) == c(7, 92)))

# Keep the raw (unstandardized) features: Problem 4 plots against raw values
TRI_raw <- TRI
TEI_raw <- TEI

# 3. Standardize each feature ROW with training mean and POPULATION sd
#    (divide by n = 300, not n - 1 like sd() does).
#    A length-7 vector is recycled down the columns of a 7 x n matrix,
#    so element i is applied to row i - exactly the per-feature transform.
TRI_mean <- rowMeans(TRI_raw)
TRI_sd   <- sqrt(rowMeans((TRI_raw - TRI_mean)^2))
TRI <- (TRI_raw - TRI_mean) / TRI_sd
TEI <- (TEI_raw - TRI_mean) / TRI_sd   # test uses TRAINING statistics

# Sanity check: every training row now has mean 0 and population sd 1
stopifnot(all(abs(rowMeans(TRI)) < 1e-10),
          all(abs(sqrt(rowMeans(TRI^2)) - 1) < 1e-10))

# Task 3.1: 8 x 8 ABSOLUTE Pearson correlation matrix (training set)
# (Pearson correlation is unchanged by standardization.)
train_df <- data.frame(t(TRI), mpg = TRO[, 1])
cor_matrix <- abs(cor(train_df))
print(round(cor_matrix, 3))

mpg_cor <- cor_matrix["mpg", feature_names]
best_feature <- names(which.max(mpg_cor))
cat("Input most correlated with mpg:", best_feature,
    sprintf("(|r| = %.3f)\n", max(mpg_cor)))

# Task 3.2: Figures (saved as PNG files for the report)

# Fig. 5 - TRO vs TRI[0,:] (cylinders, standardized)
fig5 <- ggplot(train_df, aes(x = cylinders, y = mpg)) +
  geom_point(alpha = 0.6) +
  labs(x = "Cylinders (standardized, unitless)",
       y = "Fuel economy (mpg)",
       title = "Fig. 5: Training mpg vs cylinders")
ggsave("fig5_mpg_vs_cylinders.png", fig5, width = 6, height = 4.5, dpi = 150)

# Fig. 6 - TRO vs TRI[3,:] (weight, standardized)
fig6 <- ggplot(train_df, aes(x = weight, y = mpg)) +
  geom_point(alpha = 0.6) +
  labs(x = "Weight (standardized, unitless)",
       y = "Fuel economy (mpg)",
       title = "Fig. 6: Training mpg vs weight")
ggsave("fig6_mpg_vs_weight.png", fig6, width = 6, height = 4.5, dpi = 150)

# Fig. 7 - all 7 features vs mpg in a 3 x 3 grid.
# Discrete features get a small horizontal jitter so overlapping points show.
# Jitter uses R's default Mersenne-Twister RNG with set.seed(42)
# (state this in the report; it only affects the picture, not any result).
set.seed(42)
discrete_features <- c("cylinders", "model_year", "origin")
long_df <- pivot_longer(train_df, cols = all_of(feature_names),
                        names_to = "feature", values_to = "value")
long_df$feature <- factor(long_df$feature, levels = feature_names)
is_discrete <- long_df$feature %in% discrete_features
long_df$value[is_discrete] <- long_df$value[is_discrete] +
  runif(sum(is_discrete), -0.05, 0.05)

fig7 <- ggplot(long_df, aes(x = value, y = mpg)) +
  geom_point(alpha = 0.5, size = 1) +
  facet_wrap(~ feature, ncol = 3, scales = "free_x") +
  labs(x = "Feature value (standardized, unitless)",
       y = "Fuel economy (mpg)",
       title = "Fig. 7: All seven features vs mpg (training set)")
ggsave("fig7_all_features_vs_mpg.png", fig7, width = 9, height = 8, dpi = 150)