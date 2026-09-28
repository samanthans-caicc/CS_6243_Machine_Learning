# Problem 4 (18 points): Storage Efficiency

#   1. Start with if (!exists("Kset")) source("q3.R") so you inherit Kset, r, N, and the Problem 3 error values.

if (!exists("Kset")) source("q3.R")

#   2. Compute Sfull <- 784 * N, the scalar count for the full matrix. With N = 18614 this is 14,593,376.
Sfull <- 784 * N
cat("Sfull =", Sfull, "\n")

#   3. Compute SLR <- 784 + 784 * Kset + Kset * N for every K in Kset, which simplifies to 784 + 19398 * K here.
SLR <- 784 + 19398 * Kset
cat("SLR =", SLR, "\n")

#   4. Remember the leading 784 is the mean image mu, the 784*K term is Q_K, and the K*N term is Z_K.
cat("SLR per K =", 19398, "\n")

#   5. Count scalar values only. Ignore labels, matrix dimensions, software metadata, and numerical-format
#      overhead, exactly as the PDF instructs.


rho <- SLR / Sfull
cat("rho =", rho, "\n")

#   6. Compute rho <- SLR / Sfull for every K in Kset.
rho <- SLR / Sfull
cat("rho =", rho, "\n")

#   7. Build one table with columns K, SLR(K), rho(K), and the Problem 3 error e_svd at the same K, so the
#      accuracy and storage columns sit side by side for the discussion.
table_storage <- data.frame(
  K = Kset,
  SLR = SLR,
  rho = rho,
  e_svd = e_svd
)
print(table_storage)

#   8. Note that rho(1) is about 0.0014 and rho(632) is about 0.8401, so every K in Kset saves storage.

#   9. Solve 784 + 19398*K = Sfull to find the break-even rank, which is about K = 752.
K_break_even <- (Sfull - 784) / 19398
cat("Break-even rank K =", K_break_even, "\n")
png("fig4.png")
plot(Kset, SLR, type = "l", xlab = "K", ylab = "Storage (scalars)", main = "Storage vs K")
abline(h = Sfull, lty = "dashed", col = "red")
legend("topright", legend = c("SLR", "Sfull"), col = c("black", "red"), lty = c(1, 2))
dev.off()
#  10. State that the break-even rank exceeds r = 632, so the low-rank form never costs more than the full matrix
#      anywhere in the admissible range. This is a consequence of N being much larger than 784.

#  11. Plot Fig. 4 as SLR(K) versus K and save it to its own PNG, as in q2.R.

#  12. Add a horizontal reference line at Sfull using abline(h = Sfull, lty = "dashed"), and label it.

#  13. Label the vertical axis in scalar values, since the PDF asks for that explicitly.

#  14. Since Sfull sits well above the largest SLR value, set ylim so both the curve and the reference line are
#      visible in the same frame.

#  15. In the discussion, pair each K with its Problem 3 error and point out that storage grows linearly in K
#      while the error falls fastest at small K.

#  16. Identify the useful operating range, where the error has largely flattened but rho is still small, and
#      justify it from the K90, K95, and K99 values carried over from q2.R.

#  17. Close by stating the tradeoff plainly: K buys accuracy at a constant marginal storage cost of 19398 scalars
#      per unit rank, while the accuracy return per unit rank shrinks as K grows.
