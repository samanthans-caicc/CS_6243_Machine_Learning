# -----------------------------------------------------------------------------
# PROBLEM 1 - Data splitting and standardization
# -----------------------------------------------------------------------------
# split_data(train_size)
#   - stop() with an error if train_size > 10000
#   - train = first train_size permuted indices
#   - test  = permuted positions 10001:15000 (always 5000 samples)
#   - return x_train, y_train, x_test, y_test

split_data <- function(train_size) {
  if (train_size > 10000) {
    stop("train_size must not exceed 10000")
  }

  train_indices <- perm[seq_len(train_size)]
  test_indices <- perm[10001:15000]

  list(
    x_train = x[train_indices],
    y_train = y[train_indices],
    x_test = x[test_indices],
    y_test = y[test_indices]
  )
}

# standardize(train, test)
#   - mean and POPULATION sd (divide by n, NOT sd() which divides by n-1)
#   - compute stats on train only; apply to both train and test
#   - return standardized train, standardized test, AND the stats (mu, sigma)
standardize <- function(train, test) {
  mu <- mean(train)
  sigma <- sqrt(mean((train - mu)^2))

  list(
    train = (train - mu) / sigma,
    test = (test - mu) / sigma,
    mu = mu,
    sigma = sigma
  )
}
 
# Task 1.1: print lengths of the 4 arrays for train_size = 1000, 5000, 10000.
#           Confirm test length is 5000 each time.
# Task 1.2: Fig. 0 - side-by-side scatter (train_size = 5000):
#           left = raw x vs y, right = standardized x vs y. Label axes + units.
 
