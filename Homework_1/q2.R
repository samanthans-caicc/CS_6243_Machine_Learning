#  Setup
#   1. Load Problem 0's data. Put if (!exists("X_centered")) source("q0.r") at the top of q2.R. You need X_centered
#      only; no labels this time.

if (!exists("X_centered")) source("q0.r")

#  Compute the SVD

#   2. Set q. q <- min(784, N). Since N = 18614, q = 784. Say so in the report: the economy SVD is capped by the
#      pixel dimension, not the sample count.

q <- min(784, N)  # q = 784 since N = 18614

#   3. Run the economy SVD. sv <- svd(X_centered). R's default is already economy size: sv$u is 784 x 784, sv$d is
#      length 784, sv$v is N x 784. Report those dimensions to show they match U in R^(784xq), Sigma in R^(qxq),
#      V in R^(Nxq).

sv <- svd(X_centered, nu = q, nv = q)  # economy SVD with q = 784

#   4. Keep Sigma as a vector. sv$d is the diagonal of Sigma, already sorted descending. Don't build the full
#      784 x 784 diagonal matrix; you never need it.

sigma <- sv$d  # keep Sigma as a vector of singular values
U <- sv$u      # keep the left singular vectors
V <- sv$v      # keep the right singular vectors

# Report the factor dimensions so they can be checked against U in R^(784xq), Sigma in R^(qxq), V in R^(Nxq).
cat("q =", q, "\n")
cat("dim(U):", nrow(U), "x", ncol(U), "\n")
cat("length of Sigma diagonal:", length(sigma), "(Sigma is", q, "x", q, ")\n")
cat("dim(V):", nrow(V), "x", ncol(V), "\n")

#   5. Sanity-check orthonormality (optional but cheap, and it backs up the "columns are orthonormal" claim):
#      report max(abs(crossprod(sv$u) - diag(784))).
#      Done before truncation, so it tests all q columns.

cat("max|t(U) U - I|:", format(max(abs(crossprod(U) - diag(q))), digits = 6), "\n")
cat("max|t(V) V - I|:", format(max(abs(crossprod(V) - diag(q))), digits = 6), "\n")

#  Numerical rank

#   6. Get eps_mach. eps <- .Machine$double.eps (2.220446e-16). This is the float64 machine epsilon the PDF's
#      Preprocessing section mandates.

eps <- .Machine$double.eps  # machine epsilon for float64

#   7. Compute the tolerance. tau <- max(784, N) * eps * sv$d[1], i.e. 18614 * eps_mach * sigma_1.

tau <- max(784, N) * eps * sigma[1]  # tolerance for numerical rank

#   8. Count the rank. r <- sum(sv$d > tau). Strictly greater than, per r = #{i : sigma_i > tau}.
#      Named r, not rank, both because the PDF calls it r and because `rank` shadows base::rank(). q3-q6 expect r.

r <- sum(sigma > tau)  # numerical rank of X_centered

#   9. Report tau and r explicitly; the PDF asks for both, not just r.
#      Measured result: tau = 4.863126e-07 and r = 632, well below the 784 ceiling. Two causes, in order of size:
#      (a) 124 pixel rows of X are identically zero across all 18,614 images (dead border pixels), which caps the
#          rank at 784 - 124 = 660 before anything else applies; and
#      (b) 104 rows are nonzero in fewer than 5 images, and exact dependencies among those ultra-sparse border
#          rows account for the remaining gap down to 632.
#      Mean-centering is NOT the cause; it costs at most one direction and does not bind here.

cat("Tolerance tau:", format(tau, digits = 7), "\n")
cat("Numerical rank r:", r, "\n")

#  10. Truncate. sigma <- sv$d[1:r]. Everything below uses only these r positive singular values.
#      U keeps r columns, which is exactly enough for every K in Kset, since max(Kset) = r.

sigma <- sigma[1:r]  # truncate to the positive singular values
U <- U[, 1:r]        # truncate the left singular vectors
V <- V[, 1:r]        # truncate the right singular vectors

# The full 632-value list belongs in the report figures, not the console; Q2 does not ask for it as a table.
cat("Largest 10 singular values:\n")
print(head(sigma, 10))
cat("Smallest 10 retained singular values:\n")
print(tail(sigma, 10))

# The rank is insensitive to tau: there is a ten-order-of-magnitude cliff at the cutoff, with tau inside it.
cat("sigma_r =", format(sigma[r], digits = 6),
    " sigma_(r+1) =", format(sv$d[r + 1], digits = 6), "\n")

#  Fig. 1: singular-value spectrum

#  11. Plot sigma_i vs i for i = 1..r:
#      plot(1:r, sigma, type = "l", log = "y", xlab = "Index i", ylab = expression(sigma[i]))
#      A log vertical axis is explicitly permitted and you want it; the values span several orders of magnitude.
#      Written to its own file so Fig. 1 and Fig. 2 can be captioned separately in the report.

png("fig1_singular_values.png", width = 1600, height = 1200, res = 200, type = "cairo")
plot(1:r, sigma, type = "l", log = "y", xlab = "Index i", ylab = expression(sigma[i]))
title("Singular-value spectrum of the centered data")
grid()  # add a grid for better readability
dev.off()

#  12. Caption it with the conclusion, not the contents: the spectrum decays rapidly over the first few dozen
#      components then flattens, so variance is concentrated in a low-dimensional subspace.

#  Cumulative energy

#  13. Compute E(K). energy <- cumsum(sigma^2) / sum(sigma^2), a vector of length r where energy[K] is E(K).
#      The denominator sums over the r POSITIVE singular values only, not all q.

energy <- cumsum(sigma^2) / sum(sigma^2)

#  14. Fig. 2: plot E(K) vs K. Linear axes, y from 0 to 1. Add horizontal reference lines at 0.90, 0.95, and 0.99
#      so the thresholds are readable straight off the figure.

png("fig2_cumulative_energy.png", width = 1600, height = 1200, res = 200, type = "cairo")
plot(1:r, energy, type = "l", ylim = c(0, 1), xlab = "Rank K", ylab = expression(E(K)))
title("Cumulative energy ratio of the centered data")
grid()  # add a grid for better readability
abline(h = c(0.90, 0.95, 0.99), col = "red", lty = "dashed")  # thresholds at 0.90, 0.95, 0.99
dev.off()

#  Thresholds and discussion

#  15. Find the smallest K for each threshold. K90 <- which(energy >= 0.90)[1], same for 0.95 and 0.99. Report all
#      three in a small table. Use >=, since the PDF says "at least".
#      Measured result: 61, 115, and 272.

K90 <- which(energy >= 0.90)[1]
K95 <- which(energy >= 0.95)[1]
K99 <- which(energy >= 0.99)[1]
cat("Smallest K for each threshold:\n")
print(data.frame(Threshold = c(0.90, 0.95, 0.99), K = c(K90, K95, K99)))

#  16. Discuss, covering all three terms the PDF names:
#      - Data variance: sigma_i^2 is the variance captured along direction u_i, and the total sum of sigma_i^2 is
#        ||X_centered||_F^2, so E(K) is literally the fraction of total variance retained.
#      - Redundancy: neighboring pixels are strongly correlated and 124 border pixels are identically zero, so few
#        directions carry most of the signal.
#      - Effective linear dimensionality: contrast the algebraic rank r = 632 with K90 = 61, K95 = 115, and
#        K99 = 272. The gap is the point; 61 of 632 directions carry 90% of the variance, so the data sits near a
#        subspace far smaller than its rank suggests.
#      Tie the flat tail of Fig. 1 to the plateau of Fig. 2. Explain the rank deficit using the two causes recorded
#      at step 9 (124 all-zero pixel rows, then exact dependencies among ultra-sparse border rows), and note that
#      the cutoff is unambiguous: sigma_632 = 0.934547 but sigma_633 = 7.256e-11.

#  Carry-forward

#  17. Leave sv, sigma, r, and tau in the workspace. Problems 3, 4, 5, and 6 all need Q_K <- U[, 1:K] and r for
#      Kset, so q3.R onward should source("q2.R") the way q1.R sources q0.r.
#      Note for q3: with r = 632, Kset loses 651, 701, and 751, giving {1, 51, ..., 601, 632}, i.e. 14 values.

#  Two things to watch: sv$v is 18614 x 784 (about 117 MB), so don't copy it around needlessly; and svd() on this
#  matrix takes about 30 seconds, which is another reason to compute it once here and source this file later.
