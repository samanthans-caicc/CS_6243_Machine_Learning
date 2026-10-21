# Problem 1 (4 points): Derivations
# The derivations are symbolic and belong in the report. This script only
# checks them numerically so you can trust what you write.
#
# Task 1: L(x) = (Ax - y)^T (Ax - y) = x^T A^T A x - 2 y^T A x + y^T y
#         g(x) = 2 A^T A x - 2 A^T y = 2 A^T (Ax - y) = 2(Gx - v)
# Task 2: H = d g / d x = 2 A^T A  (constant, does not depend on x)
# Task 3: g(u) - g(v) = H(u - v), so ||g(u) - g(v)|| <= ||H||_2 ||u - v||.
#         Equality holds when u - v is the top eigenvector of H, so the
#         smallest C is ||H||_2 = lambda_max(H) = 2 sigma_1^2.
# Task 4: With rank-reduced SVD A = U_r S_r V_r^T (r = rank A):
#         minimizers: x = V_r S_r^{-1} U_r^T y + z,  z in null(A) = range(V_perp)
#         L_opt = ||y||^2 - ||U_r^T y||^2 = ||(I - U_r U_r^T) y||^2
#         x_opt = V_r S_r^{-1} U_r^T y = A^+ y. Since x_opt lies in range(V_r)
#         and z is orthogonal to it, ||x||^2 = ||x_opt||^2 + ||z||^2, minimized at z = 0.
#         If A has full column rank, null(A) = {0} and the minimizer is unique.

if (!exists("x_init")) source("q0.R")

# Check Task 1: compare the analytic gradient with central finite differences
h <- 1e-6
g_fd <- sapply(1:K, function(k) {
  e <- numeric(K); e[k] <- h
  (L_fun(x_init + e) - L_fun(x_init - e)) / (2 * h)
})
cat("Max |analytic - finite-difference gradient| at x_init:",
    format(max(abs(g_fun(x_init) - g_fd)), digits = 4), "\n")

# Check Task 3: C = ||H||_2 equals 2 * sigma_1^2
cat("lambda_max(2 A^T A) vs 2 sigma_1^2:",
    max(eigen(2 * G, symmetric = TRUE, only.values = TRUE)$values),
    2 * svd(A)$d[1]^2, "\n")
