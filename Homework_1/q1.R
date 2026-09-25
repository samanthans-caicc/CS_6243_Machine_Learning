#  Setup
#   1. Load Problem 0's data. Put if (!exists("X_centered")) source("q0.r") at the top of q1.R. You'll use X_centered and
#      labels.

if (!exists("X_centered")) source("q0.r")

#   Compute

#   2. Loop over the classes. Go through c = 0, 1, 9.
classes <- c(0, 1, 9)
centroids <- list()
norms <- list()

#   3. Pick that class's columns. X̃_c is X_centered restricted to the columns where labels == c.
for (c in classes) {
  X_c <- X_centered[, labels == c]
  N_c <- ncol(X_c)
  m_c <- rowSums(X_c) / N_c
  centroids[[as.character(c)]] <- m_c
  
  l1 <- sum(abs(m_c))
  l2 <- sqrt(sum(m_c^2))
  l3 <- (sum(abs(m_c)^3))^(1/3)
  norms[[as.character(c)]] <- c(l1 = l1, l2 = l2, l3 = l3)
}

#   4. Count them. N_c is the number of columns in X̃_c. You should get 5923, 6742, and 5949.
cat("Number of images per class:\n")
for (c in classes) {
  cat("Class", c, ":", sum(labels == c), "\n")
}

#   5. Compute the centroid. m_c = (1/N_c)·X̃_c·1, which is the average of each pixel (each row) across that class's images.
#      It has 784 values, the same kind of calculation as μ in q0.
cat("Centroids per class:\n")
for (c in classes) {
  cat("Class", c, "centroid (first 10 values):", head(centroids[[as.character(c)]], 10), "\n")
}

#   6. Compute three norms for each m_c:
#      - ℓ1: the sum of the absolute values
#      - ℓ2: the square root of the sum of squares
#      - ℓ3: the cube root of the sum of |values|³

#      Use the p-norm formula from page 1 of the PDF directly, since base R has no built-in ℓ3 norm.

norm_df <- as.data.frame(do.call(rbind, norms))
#   Present

#   7. Build one labeled table. Use one row per class (0, 1, 9) and one column per norm. A data.frame works; print it with
#      enough digits, e.g. print(df, digits = 10).

cat("Labeled table of norms:\n")
print(norm_df, digits = 10)

