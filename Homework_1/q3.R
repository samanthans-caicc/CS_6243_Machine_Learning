# Problem 3 (18 points): Low-Rank Approximation

#   1. Start with if (!exists("sigma", inherits = FALSE)) source("q2.R") so you inherit U, sigma, r, X_centered, mu, and N.

# inherits = FALSE matters: stats::sigma() is a real function, so a bare
# exists("sigma") is always TRUE and this guard would never fire.
if (!exists("sigma", inherits = FALSE)) source("q2.R")

#   2. Build Kset <- sort(unique(c(intersect(seq(1, 751, by = 50), 1:r), r))). With r = 632 this gives 14 values:
#      1, 51, 101, 151, 201, 251, 301, 351, 401, 451, 501, 551, 601, 632.
Kset <- sort(unique(c(intersect(seq(1, 751, by = 50), 1:r), r)))

#   3. Print Kset and state in the report that 651, 701, and 751 were dropped because they exceed r = 632.\
cat("Kset =", Kset, "\n")

#   4. Compute the SVD-based error for every K with tail <- rev(cumsum(rev(sigma^2))), then
#      e_svd <- sqrt(c(tail[-1], 0))[Kset]. Entry K is the square root of the discarded tail above K.
tail <- rev(cumsum(rev(sigma^2)))
e_svd <- sqrt(c(tail[-1], 0))[Kset]

#   5. Check that e_svd at K = r is exactly 0, since the empty sum is zero by the convention stated in the PDF.
cat("e_svd at K = r =", e_svd[which(Kset == r)], "\n")

#   6. Also compute the SVD-based error at K = 100 separately. The PDF requires it for verification even though
#      100 is not in Kset (Kset contains 101, not 100).
e_svd_100 <- sqrt(c(tail[-1], 0))[100]
cat("e_svd at K = 100 =", e_svd_100, "\n")

e_svd_1 <- sqrt(c(tail[-1], 0))[1]
cat("e_svd at K = 1 =", e_svd_1, "\n")

#   7. Write a small function of K that forms Q <- U[, 1:K, drop = FALSE] and Y <- Q %*% (crossprod(Q, X_centered)).
svd_error <- function(K) {
  Q <- U[, 1:K, drop = FALSE]
  Y <- Q %*% (crossprod(Q, X_centered))
  return(sqrt(sum((X_centered - Y)^2)))
}

#   8. Keep the parentheses in that order so you multiply 784 x K by K x N and never form anything larger than X.
direct_error <- function(K) {
  Q <- U[, 1:K, drop = FALSE]
  Y <- Q %*% (crossprod(Q, X_centered))
  return(sqrt(sum((X_centered - Y)^2)))
}

#   9. Compute the direct error norm(X_centered - Y, "F") at K = 1, K = 100, and K = r only, since those are the
#      three values the PDF asks you to verify directly.
direct_error_1 <- direct_error(1)
cat("direct error at K = 1 =", direct_error_1, "\n")
direct_error_100 <- direct_error(100)
cat("direct error at K = 100 =", direct_error_100, "\n")
direct_error_r <- direct_error(r)
cat("direct error at K = r =", direct_error_r, "\n")

#  10. Do not store all the Y matrices. Each is about 117 MB, so compute the error inside the loop and let Y go.

#  11. Form X_hat <- mu %*% t(rep(1, N)) + Y for at least one K so you have built it, and say in the report that
#      X_hat is the affine reconstruction in original image coordinates while Y is the linear one in centered
#      coordinates.
K_example <- 100  # choose an example K for which to form X_hat
Q_example <- U[, 1:K_example, drop = FALSE]
Y_example <- Q_example %*% (crossprod(Q_example, X_centered))
X_hat <- mu %*% t(rep(1, N)) + Y_example

#  12. Compare direct against SVD-based at K = 1, 100, and r, reporting the absolute and relative difference.
abs_diff_1 <- abs(direct_error_1 - e_svd_1)
rel_diff_1 <- abs_diff_1 / direct_error_1
cat("Absolute difference at K = 1:", abs_diff_1, "\n")
cat("Relative difference at K = 1:", rel_diff_1, "\n")

abs_diff_100 <- abs(direct_error_100 - e_svd_100)
rel_diff_100 <- abs_diff_100 / direct_error_100
cat("Absolute difference at K = 100:", abs_diff_100, "\n")
cat("Relative difference at K = 100:", rel_diff_100, "\n")

abs_diff_r <- abs(direct_error_r - e_svd[which(Kset == r)])
rel_diff_r <- abs_diff_r / direct_error_r
cat("Absolute difference at K = r:", abs_diff_r, "\n")
cat("Relative difference at K = r:", rel_diff_r, "\n")

#  13. Expect agreement to roughly 1e-10 relative rather than exact equality, and say so rather than claiming the
#      values are identical.


#  14. Plot Fig. 3 as e_svd versus Kset and save it to its own PNG, the way q2.R writes its two figures.
png("fig3.png")
plot(Kset, e_svd, type = "l", xlab = "K", ylab = "Error", main = "SVD Error vs K")
points(c(1, 100, r), c(direct_error_1, direct_error_100, direct_error_r), col = "red", pch = 19)
dev.off()

#  15. Overlay the three direct values as points on Fig. 3. The PDF allows this and it makes the verification
#      visible in the figure itself.
png("fig3_overlay.png")
plot(Kset, e_svd, type = "l", xlab = "K", ylab = "Error", main = "SVD Error vs K with Direct Errors")
points(c(1, 100, r), c(direct_error_1, direct_error_100, direct_error_r), col = "red", pch = 19)
dev.off()

#  16. In the discussion, say the error at rank K is the square root of the discarded tail of sigma^2, so the curve
#      falls fast where the spectrum falls fast and flattens once only tiny singular values remain.

#  17. Close the discussion with Eckart-Young: no rank-K matrix beats Q_K Q_K^T X in Frobenius norm, and the
#      minimum achievable error is exactly that discarded tail.
