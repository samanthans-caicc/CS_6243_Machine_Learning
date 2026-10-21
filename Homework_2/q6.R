# Problem 6 (15 points): Comparison of the Three Methods
# FSS p = 1, ELS, BLS (c = 0.1, beta = 0.5). Runtimes exclude CSV loading,
# G/v/C precomputation and the SVD benchmark (those happen in q0.R / q2.R).

if (!exists("els_run")) source("q4.R")
if (!exists("table1"))  source("q5.R")

cmp <- list(
  "FSS (p = 1)"             = function() gd_run("fss", gamma = 1 / C),
  "ELS"                     = function() gd_run("els"),
  "BLS (c = 0.1, beta = 0.5)" = function() gd_run("bls", beta = 0.5, c_arm = 0.1)
)
cmp_runs  <- lapply(cmp, function(f) f())
cmp_times <- lapply(cmp, time_runs)
cmp_gap   <- sapply(cmp_runs, obj_path) - L_opt     # 30 x 3, raw values

# ---- Task 1: Fig. 6 ----
# Log scale: gaps that are zero or slightly negative from roundoff cannot be
# drawn, so ONLY for plotting they are clipped to a floor of
# eps * L_opt (the roundoff level of L). Tables use the raw values.
gap_floor <- .Machine$double.eps * L_opt
gap_plot  <- pmax(cmp_gap, gap_floor)
cat("\nNumber of gap values clipped for plotting:", sum(cmp_gap < gap_floor), "\n")

cmp_cols <- c("#d95f02", "blue", "#7570b3")
png("fig6_comparison.png", width = 1600, height = 1100, res = 200)
matplot(1:n_iter, gap_plot, type = "l", lty = 1, lwd = 2, col = cmp_cols, log = "y",
        xlab = "Iteration n", ylab = expression(L(x[n]) - L[opt] ~ "(log scale)"),
        main = "Fig. 6: Objective gap of the three GD methods")
abline(h = gap_floor, lty = 3, col = "grey40")
legend("topright", legend = c(names(cmp), "plot floor (eps * L_opt)"),
       col = c(cmp_cols, "grey40"), lty = c(1, 1, 1, 3), lwd = 2, bty = "n")
dev.off()

# ---- Task 2: Table 2 ----
table2 <- data.frame(
  method     = names(cmp),
  final_gap  = cmp_gap[n_iter, ],
  updates    = sapply(cmp_runs, `[[`, "n_updates"),
  obj_evals  = sapply(cmp_runs, `[[`, "n_evals"),
  median_us  = sapply(cmp_times, median) * 1e6,
  iqr_us     = sapply(cmp_times, IQR) * 1e6,
  reps       = sapply(cmp_times, length),
  row.names  = NULL
)
cat("\nTable 2 (runtimes in microseconds; exclude preprocessing and SVD):\n")
print(table2, digits = 5)

# Work per iteration NOT in the objective-evaluation count (K = 50):
#   gradient 2(Gx - v):        one K x K matvec       (all methods)
#   ELS step ||g||^2/(2 g^T G g): one extra K x K matvec (ELS only)
#   each BLS objective eval:   one N x K matvec        (BLS only)
cat("\nPer-iteration matvecs: FSS = 1 (KxK); ELS = 2 (KxK);",
    "BLS = 1 (KxK) + (1 + backtracks) (NxK)\n")

# Gap reached after each method's n-th iteration, handy for 'progress per iteration'
cat("\nGap at n = 1, 5, 10, 20, 30:\n")
print(signif(cmp_gap[c(1, 5, 10, 20, 30), ], 4))

# Save all timings for the variability discussion
all_times <- do.call(rbind, lapply(names(cmp_times), function(m)
  data.frame(method = m, rep = seq_along(cmp_times[[m]]), seconds = cmp_times[[m]])))
write.csv(all_times, "timings_problem6.csv", row.names = FALSE)
