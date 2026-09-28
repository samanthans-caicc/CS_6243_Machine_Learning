# Problem 6 (18 points): Low-Rank Denoising

#   1. Start with if (!exists("Kset")) source("q3.R") so you inherit Kset, U, mu, X, X_centered, N, and r.

if (!exists("Kset")) source("q3.R")

#   2. Use the original uncentered X from Problem 0 as the clean matrix. Do not use X_centered as the starting
#      point for the corruption.

#   3. Call set.seed(4463) once, immediately before generating any random values.
set.seed(4463)
# (reported below, after the draws)

#   4. Report the generator with RNGkind(), which on stock R 4.3.3 returns Mersenne-Twister, Inversion, Rejection.
#      The PDF requires naming the algorithm because realizations differ across environments.

#   5. Draw the corruption mask as M <- matrix(runif(784 * N) < 0.05, nrow = 784), which is 14.6 million draws.

M <- matrix(runif(784 * N) < 0.05, nrow = 784)

#   6. Draw the replacement values as B <- matrix(runif(784 * N) < 0.5, nrow = 784) in a second independent call.

B <- matrix(runif(784 * N) < 0.5, nrow = 784)

#   7. Build Xn <- X, then set Xn[M] <- ifelse(B[M], 255, 0) so only masked entries change.

Xn <- X
Xn[M] <- ifelse(B[M], 255, 0)

#   8. Leave every unmasked pixel exactly as it was. The PDF forbids clipping or otherwise modifying uncorrupted
#      pixels.

#   9. Verify the corruption rate by reporting mean(M), which should sit near 0.05.
#      Do NOT expect Xn to differ from X in 5 percent of entries: most MNIST pixels are already 0, so a masked
#      background pixel that draws B = 0 is overwritten with the value it already had. The observed differ-rate
#      is about 0.029, and that is correct behaviour, not a bug.

cat("mean(M), target 0.05:", mean(M), "\n")
cat("fraction of entries actually changed:", sum(Xn != X) / (784 * N), "\n")
cat("  (below 0.05 because masked background pixels often redraw to 0)\n")
cat("RNG algorithm:", RNGkind(), "\n")

#  10. Center the noisy matrix with the clean mean from Problem 0: Xn_centered <- sweep(Xn, 1, mu). Do not
#      recompute a mean from the noisy data.

Xn_centered <- sweep(Xn, 1, mu)

#  11. For every K in Kset, project using Q_K from the clean data: Xn_hat <- mu %*% t(rep(1, N)) +
#      Q %*% crossprod(Q, Xn_centered), with Q <- U[, 1:K, drop = FALSE].

#  14. Compute the constant e_baseline <- norm(X - Xn, "F") once, outside the loop.

e_baseline <- norm(X - Xn, "F")

#  12. Keep the parentheses so the K x N product forms first and nothing larger than X is built.

#  13. Compute e_clean <- norm(X - Xn_hat, "F") and e_noisy <- norm(Xn - Xn_hat, "F") inside the loop, then
#      discard Xn_hat before the next K.
#  15. Do not clip the denoised reconstructions before computing any of these errors.

# One pass per K: build the reconstruction, take both errors, then drop it.
# Xn_hat is 784 x 18614 (about 117 MB), so it must not be accumulated.
mu1 <- mu %*% t(rep(1, N))
e_clean <- numeric(length(Kset))
e_noisy <- numeric(length(Kset))

for (i in seq_along(Kset)) {
  K <- Kset[i]
  Q <- U[, 1:K, drop = FALSE]          # clean-data subspace, from q2.R
  Xn_hat <- mu1 + Q %*% crossprod(Q, Xn_centered)
  e_clean[i] <- norm(X  - Xn_hat, "F")
  e_noisy[i] <- norm(Xn - Xn_hat, "F")
  cat("K =", K, " e_clean =", e_clean[i], " e_noisy =", e_noisy[i], "\n")
  rm(Xn_hat)
}
rm(mu1)

#  16. Plot Fig. 7 with all three errors against K on one axis, with e_baseline drawn as a horizontal reference
#      line, and save it to its own PNG.
#  17. Add a legend naming the three curves, since the figure is unreadable without one.

png("fig7_denoising.png", width = 1600, height = 1200, res = 200, type = "cairo")
plot(Kset, e_clean, type = "b", pch = 19, col = "blue",
     ylim = range(0, e_clean, e_noisy, e_baseline),
     ylab = "Frobenius error", xlab = "Rank K",
     main = "Low-rank denoising error versus rank")
lines(Kset, e_noisy, type = "b", pch = 17, col = "red")
# e_baseline is a single number, so it is a reference LINE, not a curve.
abline(h = e_baseline, col = "darkgreen", lty = "dashed", lwd = 2)
grid()
legend("right",
       legend = c("e_clean(K)", "e_noisy(K)", "e_baseline"),
       col = c("blue", "red", "darkgreen"),
       pch = c(19, 17, NA), lty = c(1, 1, 2), bty = "n")
dev.off()

#  18. Report the K that minimizes e_clean, and also report a near-minimizing range rather than only the single
#      argmin, because the curve is typically flat near its bottom.

K_min_e_clean <- Kset[which.min(e_clean)]
K_min_range   <- Kset[e_clean <= min(e_clean) * 1.05]

cat("\ne_baseline:", e_baseline, "\n")
cat("K minimizing e_clean:", K_min_e_clean, "\n")
cat("min e_clean:", min(e_clean), "\n")
cat("K within 5% of the minimum:", K_min_range, "\n")
cat("fraction of baseline corruption remaining:", min(e_clean) / e_baseline, "\n")

print(data.frame(K = Kset, e_clean = e_clean, e_noisy = e_noisy), row.names = FALSE)

# -------------------------------
#           DISCUSSION
# --------------------------------

#  19. In the discussion, explain that e_noisy falls monotonically in K because a larger subspace fits the noisy
#      matrix better, so it is not a denoising criterion.

#  20. Explain that e_clean falls and then rises, giving a U shape, because early directions restore digit
#      structure while later directions start fitting the salt-and-pepper corruption.

#  21. Explain that a small rank discards meaningful digit structure along with the noise, which is why e_clean is
#      high at K = 1 despite the subspace containing almost no corruption.

#  22. Compare the minimum of e_clean against e_baseline to state how much of the corruption the projection
#      actually removed, and say so as a ratio.

#  23. Connect the minimizer back to the K90, K95, and K99 values from q2.R and note whether denoising favors a
#      smaller rank than pure reconstruction does.
