# -----------------------------------------------------------------------------
# PROBLEM 5 - Multivariate polynomial regression, all 7 features
# -----------------------------------------------------------------------------
# Requires, in order: setup.R, q1.R, q2.R (standardize_poly, solve_LS, mse),
# q3.R (standardized TRI, TEI and TRO, TEO), q4.R (test_mse_mat, for 5.4).

library(ggplot2)

# 1. All exponent vectors (e1, ..., e7) with e_j >= 0 and sum(e) <= M.
#    Rows are sorted by total degree, then so that x1 terms come first.
#    Row 1 is always all zeros = the bias term.
exponent_vectors <- function(M, p = 7) {
  E <- as.matrix(expand.grid(rep(list(0:M), p)))
  E <- E[rowSums(E) <= M, , drop = FALSE]
  ord <- do.call(order, c(list(rowSums(E)), lapply(seq_len(p), function(j) -E[, j])))
  E <- E[ord, , drop = FALSE]
  dimnames(E) <- NULL
  E
}

# Readable name for one exponent vector, e.g. c(2,1,0,0,0,0,0) -> "x1^2*x2"
term_name <- function(e) {
  if (all(e == 0)) return("1")
  parts <- ifelse(e == 1, paste0("x", seq_along(e)),
                  paste0("x", seq_along(e), "^", e))
  paste(parts[e > 0], collapse = "*")
}

# 2. Design matrix (n x number_of_terms). Z is n x 7 (samples in rows).
#    Column k = product over j of Z[, j]^E[k, j].
multi_poly_features <- function(Z, E) {
  X <- matrix(1, nrow = nrow(Z), ncol = nrow(E))
  for (k in seq_len(nrow(E))) {
    for (j in which(E[k, ] > 0)) {
      X[, k] <- X[, k] * Z[, j]^E[k, j]
    }
  }
  X
}

# Inputs: the STANDARDIZED features from q3.R, transposed to samples-in-rows
Z_train <- t(TRI)   # 300 x 7
Z_test  <- t(TEI)   #  92 x 7

degrees_multi <- 1:3
results <- data.frame(M = degrees_multi, n_terms = NA_integer_,
                      train_mse = NA_real_, test_mse = NA_real_)

for (i in seq_along(degrees_multi)) {
  M <- degrees_multi[i]
  E <- exponent_vectors(M)
  stopifnot(nrow(E) == choose(7 + M, M))   # 8, 36, 120

  # Task 5.1: explicit term list
  terms <- apply(E, 1, term_name)
  cat(sprintf("\nM = %d: %d terms\n", M, length(terms)))
  cat(paste(terms, collapse = ", "), "\n")

  # 3. Build, column-standardize (bias skipped), solve, evaluate
  X_train <- multi_poly_features(Z_train, E)
  X_test  <- multi_poly_features(Z_test, E)
  ps <- standardize_poly(X_train, X_test)
  w <- solve_LS(ps$X_train, TRO[, 1])

  results$n_terms[i]   <- length(terms)
  results$train_mse[i] <- mse(ps$X_train %*% w, TRO[, 1])
  results$test_mse[i]  <- mse(ps$X_test %*% w, TEO[, 1])
}

# Task 5.2: table of train/test MSE (mpg^2)
cat("\n"); print(transform(results, train_mse = round(train_mse, 3),
                                    test_mse  = round(test_mse, 3)))

# Task 5.3: Fig. 11 - training (red), testing (blue)
mse_long <- rbind(data.frame(M = results$M, mse = results$train_mse, set = "Training"),
                  data.frame(M = results$M, mse = results$test_mse,  set = "Testing"))
fig11 <- ggplot(mse_long, aes(x = M, y = mse, colour = set)) +
  geom_line() +
  geom_point(size = 2) +
  scale_colour_manual(values = c(Training = "red", Testing = "blue")) +
  scale_x_continuous(breaks = degrees_multi) +
  labs(x = "Total polynomial degree M", y = "MSE (mpg^2)", colour = NULL,
       title = "Fig. 11: Multivariate polynomial regression, MSE vs M")
ggsave("fig11_multivariate_mse_vs_M.png", fig11, width = 6.5, height = 4.5, dpi = 150)

# Task 5.4: numbers for the discussion (the discussion itself goes in the report)
best_multi <- results[which.min(results$test_mse), ]
cat(sprintf("\nBest multivariate: M = %d (%d terms), test MSE = %.3f\n",
            best_multi$M, best_multi$n_terms, best_multi$test_mse))
if (exists("test_mse_mat")) {
  best_uni <- which(test_mse_mat == min(test_mse_mat), arr.ind = TRUE)[1, ]
  cat(sprintf("Best univariate (Problem 4): %s, M = %d, test MSE = %.3f\n",
              rownames(test_mse_mat)[best_uni[1]], best_uni[2] - 1, min(test_mse_mat)))
}