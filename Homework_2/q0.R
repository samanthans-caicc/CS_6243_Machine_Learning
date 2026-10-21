# Problem 0 (0 points): Setup
# Run every script from inside Homework_2/ (setwd there first). Each later
# script sources the previous one if its objects are missing, like Homework 1.

# Load A, y, x_init. The CSVs are whitespace separated with no header.
A      <- as.matrix(read.table("A.csv", header = FALSE))
y      <- as.numeric(read.table("y.csv", header = FALSE)[, 1])
x_init <- as.numeric(read.table("x_init.csv", header = FALSE)[, 1])
dimnames(A) <- NULL

N <- nrow(A)   # number of rows    (100)
K <- ncol(A)   # number of columns (50)

# R stores every numeric as an IEEE double, so double precision is automatic
stopifnot(is.double(A), is.double(y), is.double(x_init))
stopifnot(length(y) == N, length(x_init) == K)

# Common precomputation used by all GD methods (lecture implementation):
# G = A^T A, v = A^T y, and g(x) = 2(Gx - v)
G <- crossprod(A)        # K x K
v <- drop(crossprod(A, y))

# Objective via the residual expression L(x) = ||Ax - y||^2
L_fun <- function(x) sum((A %*% x - y)^2)

# Gradient via the precomputed quantities
g_fun <- function(x) 2 * drop(G %*% x - v)

# One gradient-norm stopping tolerance used by every method (Problems 3-6).
# It is checked at x_{n-1} before computing step n.
grad_tol <- 1e-8

# Number of GD updates per run
n_iter <- 30

# Reproducibility info for the report
cat("N =", N, " K =", K, "\n")
cat("R version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")
cat("BLAS:", extSoftVersion()[["BLAS"]], "\n")
cat("Gradient-norm stopping tolerance:", grad_tol, "\n")
cat("No randomness is used, so no seed is needed.\n")
