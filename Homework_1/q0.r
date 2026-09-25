# Problem 0 tasks

#  Load the data. Load the MNIST training images (60,000) and their labels.

source("mnist.R")
data <- load_mnist()
images <- data$images
labels <- data$labels

# Filter images and labels to keep only 0, 1, and 9
keep <- labels %in% c(0, 1, 9)
images <- images[, , keep]
labels <- labels[keep]
original_indices <- which(keep) - 1  # 0-based indices
N <- length(labels)

# Convert to float (R only has float64). storage.mode keeps the 28 x 28 x N
# dimensions; as.numeric() would flatten the array into one long vector.
storage.mode(images) <- "double"

# Vectorize: flatten each 28x28 image row by row into a 784-long column
X <- matrix(0, nrow = 784, ncol = N)
for (i in 1:N) {
  # images[, , i] is stored transposed (see mnist.R), so column-major
  # flattening here already gives Python's row-major image.reshape(-1)
  X[, i] <- as.numeric(images[, , i])
}

# Keep labels aligned
y <- labels

# Compute the mean image: mu = (1/N) X 1, one average per pixel (784 values)
mu <- rowSums(X) / N

# Center the data: subtract mu from every column (margin 1 = rows)
X_centered <- sweep(X, 1, mu)

# Verify: X_centered %*% 1 sums each row, so use rowSums
center_check <- rowSums(X_centered)
max_center_check <- max(abs(center_check))

# Report the results
cat("N:", N, "\n")
cat("Dimensions of X:", nrow(X), "x", ncol(X), "\n")
cat("Max absolute entry of X_centered %*% 1:", format(max_center_check, digits = 6), "\n")


#   2. Filter. Keep only images labeled 0, 1, or 9, in their original order.
#   3. Save the original indices. Record each kept image's 0-based position in the 60,000-image training set. You need these for Problem 5.
#   4. Count. Set N to the number of images kept. It should be 18,614.
#   5. Convert to float. Change the pixels to float32, then to float64. In R this is just as.numeric(), since R only has float64. Note that in
#      the report.
#   6. Vectorize. Flatten each 28×28 image row by row into a 784-long column.
#   7. Build X. Stack the columns into X, a 784 × N matrix.
#   8. Keep labels aligned. Keep a label vector so column i of X still matches its digit.
#   9. Compute the mean image. μ = (1/N)·X·1, a 784-long vector.
#   10. Center the data. X̃ = X − μ1ᵀ, which subtracts μ from every column.
#   11. Verify. Compute X̃·1 and find its largest absolute entry. It should be close to zero, around 1e-10 or smaller.
#   12. Report. Report N, the dimensions of X (784 × 18614), and the max absolute entry from task 11.