# Problem 5 (30 points): GD with Backtracking Line Search
# Trial step resets to gamma = 1 each iteration; accept the first gamma with
#   L(x - gamma g) <= L(x) - c gamma ||g||^2   (equality accepted)

if (!exists("gd_run")) source("q3.R")

bls_cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")

# Two-panel figure helper: (a) objective, (b) accepted steps
plot_bls_panels <- function(file, runs, labels, cols, title) {
  obj <- sapply(runs, obj_path)
  gam <- sapply(runs, `[[`, "gam")
  png(file, width = 2400, height = 1000, res = 200)
  par(mfrow = c(1, 2))
  matplot(1:n_iter, obj, type = "l", lty = 1, lwd = 2, col = cols, log = "y",
          xlab = "Iteration n", ylab = expression(L(x[n]) ~ "(log scale)"),
          main = paste0("(a) ", title, ": objective"),
          ylim = range(c(obj, L_opt)))
  abline(h = L_opt, lty = 2, lwd = 2)
  legend("topright", legend = c(labels, "L_opt"), col = c(cols, "black"),
         lty = c(rep(1, length(cols)), 2), lwd = 2, bty = "n")
  matplot(1:n_iter, gam, type = "o", pch = 20, lty = 1, lwd = 1.5, col = cols, log = "y",
          xlab = "Iteration n", ylab = expression(gamma[n] ~ "(log scale)"),
          main = paste0("(b) ", title, ": accepted step sizes"),
          ylim = range(c(gam, 1 / C), na.rm = TRUE) * c(1, 5))  # headroom for legend
  abline(h = 1 / C, lty = 2, lwd = 2)
  legend("top", legend = c(labels, "1/C"), col = c(cols, "black"), ncol = 3,
         lty = c(rep(1, length(cols)), 2), lwd = 2, bty = "n")
  dev.off()
}

# ---- Task 1: effect of beta (c = 0.1) ----
beta_set <- c(0.1, 0.2, 0.5, 0.9)
bls_beta <- lapply(beta_set, function(b) gd_run("bls", beta = b, c_arm = 0.1))
names(bls_beta) <- paste0("beta=", beta_set)
plot_bls_panels("fig4_bls_beta.png", bls_beta, paste("beta =", beta_set),
                bls_cols, "BLS, c = 0.1")

cat("\nAccepted steps x C (first 10) for each beta:\n")
for (nm in names(bls_beta)) cat(sprintf("%-9s", nm), round(bls_beta[[nm]]$gam[1:10] * C, 3), "\n")
# Accepted steps can only be beta^k. With gamma starting at 1 >> 1/C, the
# accepted step is the first power of beta below the Armijo threshold, so
# small beta overshoots downward (tiny steps), large beta lands close to it.
cat("Powers of beta near 1/C =", 1 / C, ":\n")
for (b in beta_set) {
  k <- ceiling(log(1 / C) / log(b))
  cat(sprintf("beta=%.1f: beta^%d = %.3e, beta^%d = %.3e\n", b, k - 1, b^(k - 1), k, b^k))
}

# ---- Task 2: effect of c (beta = 0.5) ----
c_set <- c(0.01, 0.1)
bls_c <- lapply(c_set, function(cc) gd_run("bls", beta = 0.5, c_arm = cc))
names(bls_c) <- paste0("c=", c_set)
plot_bls_panels("fig5_bls_c.png", bls_c, paste("c =", c_set),
                bls_cols[c(1, 3)], "BLS, beta = 0.5")

# Report note: with this data the two settings are NOT identical (final gaps
# ~6.2e-4 for c = 0.01 vs ~1.1e-4 for c = 0.1), so discuss the difference
# rather than the "identical results" case.
cat("\nIdentical accepted steps for c = 0.01 and c = 0.1?",
    isTRUE(all.equal(bls_c[[1]]$gam, bls_c[[2]]$gam)), "\n")
# Why they can match: for a quadratic, the trial is accepted iff
# gamma <= 2(1 - c) / RQ_H(g), with RQ_H(g) = g^T H g / g^T g. Changing c from
# 0.1 to 0.01 only moves this threshold by the factor 0.99/0.9 = 1.1; if no
# power of 0.5 falls between the two thresholds, the accepted step is identical.
thr <- function(run, cc) {
  sapply(1:run$n_updates, function(n) {
    g <- g_fun(run$X[, n]); 2 * (1 - cc) * sum(g^2) / sum(g * (H %*% g))
  })
}
cat("Armijo thresholds (first 5), c=0.01:", signif(thr(bls_c[[1]], 0.01)[1:5], 4), "\n")
cat("Armijo thresholds (first 5), c=0.1 :", signif(thr(bls_c[[2]], 0.1)[1:5], 4), "\n")

# ---- Task 3: Table 1 (distinct settings; (0.1, 0.5) appears once) ----
# Report note: runtimes vary slightly between runs, so state which run the
# reported medians came from (the IQR column shows the variability).
settings <- unique(rbind(data.frame(c = 0.1, beta = beta_set),
                         data.frame(c = c_set, beta = 0.5)))
table1 <- do.call(rbind, lapply(seq_len(nrow(settings)), function(i) {
  cc <- settings$c[i]; bb <- settings$beta[i]
  run   <- gd_run("bls", beta = bb, c_arm = cc)
  times <- time_runs(function() gd_run("bls", beta = bb, c_arm = cc))
  data.frame(c = cc, beta = bb,
             final_gap   = tail(obj_path(run), 1) - L_opt,
             updates     = run$n_updates,
             obj_evals   = run$n_evals,
             backtracks  = run$n_backtracks,
             median_us   = median(times) * 1e6,
             iqr_us      = IQR(times) * 1e6,
             reps        = length(times))
}))
cat("\nTable 1 (runtimes in microseconds, ", n_rep, " timed reps after 1 warm-up):\n", sep = "")
print(table1, digits = 5)
