# Problem 2 (6 points): Numerical Computations

if (!exists("x_init")) source("q0.R")

# ---- Task 1: H and C ----
H <- 2 * G
eig_H <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
C <- max(eig_H)   # smallest global gradient-Lipschitz constant, ||H||_2

cat("H is", nrow(H), "x", ncol(H), "\n")
cat("H[1:5, 1:5] (put the full matrix in an appendix or describe it):\n")
print(round(H[1:5, 1:5], 4))
cat("C = lambda_max(H) =", format(C, digits = 10), "\n")
write.csv(H, "H.csv", row.names = FALSE)

# ---- Task 2: x_opt and L_opt from the SVD ----
sv <- svd(A)              # A = U diag(d) V^T (thin SVD)
d  <- sv$d

# Numerical rank tolerance, same rule as numpy.linalg.matrix_rank:
# tol = max(N, K) * eps * sigma_1
rank_tol <- max(N, K) * .Machine$double.eps * d[1]
r <- sum(d > rank_tol)

U_r <- sv$u[, 1:r, drop = FALSE]
V_r <- sv$v[, 1:r, drop = FALSE]
d_r <- d[1:r]

# x_opt = V_r S_r^{-1} U_r^T y
Uty   <- drop(crossprod(U_r, y))
x_opt <- drop(V_r %*% (Uty / d_r))

# Two equivalent ways to get L_opt; they should agree to roundoff
L_opt         <- L_fun(x_opt)
L_opt_formula <- sum(y^2) - sum(Uty^2)

cat("Singular values: max =", d[1], " min =", d[K], "\n")
cat("Condition number sigma_1 / sigma_K =", d[1] / d[K], "\n")
cat("Rank tolerance =", format(rank_tol, digits = 4), "  numerical rank r =", r, "\n")
cat("Rank deficient?", r < K, "\n")
cat("L_opt (residual)  =", format(L_opt, digits = 15), "\n")
cat("L_opt (formula)   =", format(L_opt_formula, digits = 15), "\n")
cat("||x_opt||_2 =", sqrt(sum(x_opt^2)), "\n")
cat("x_opt =\n"); print(round(x_opt, 6))

# Positive eigenvalues of H (used in Problem 3 Task 3). Eigenvalues of H are 2 d^2.
lam_pos <- 2 * d_r^2
lam_min_pos <- min(lam_pos)
cat("Smallest positive eigenvalue of H =", lam_min_pos, "\n")
