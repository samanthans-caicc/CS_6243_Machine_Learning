# Problem 4 (20 points): GD with Exact Line Search
# Closed-form LS step: minimize phi(gamma) = ||A(x - gamma g) - y||^2
#   phi'(gamma) = -g^T g + 2 gamma g^T G g = 0  =>  gamma = ||g||^2 / (2 g^T G g)

if (!exists("gd_run")) source("q3.R")

els_run <- gd_run("els")
els_obj <- obj_path(els_run)

cat("\nELS status:", els_run$status, "\n")
cat("ELS updates performed:", els_run$n_updates, "\n")
cat("ELS final gap L(x_30) - L_opt =", format(els_obj[n_iter] - L_opt, digits = 6), "\n")

# ---- Task 1: Fig. 2 ----
png("fig2_els_objective.png", width = 1600, height = 1100, res = 200)
plot(1:n_iter, els_obj, type = "l", lwd = 2, col = "blue", log = "y",
     xlab = "Iteration n", ylab = expression(L(x[n]) ~ "(log scale)"),
     main = "Fig. 2: GD with exact line search",
     ylim = range(c(els_obj, L_opt)))
abline(h = L_opt, lty = 2, lwd = 2, col = "black")
legend("topright", legend = c("ELS", expression(L[opt])),
       col = c("blue", "black"), lty = c(1, 2), lwd = 2, bty = "n")
dev.off()

# ---- Task 3: Fig. 3 ----
png("fig3_els_steps.png", width = 1600, height = 1100, res = 200)
plot(1:n_iter, els_run$gam, type = "o", pch = 20, lwd = 2, col = "blue",
     xlab = "Iteration n", ylab = expression(gamma[n]),
     main = "Fig. 3: ELS accepted step sizes",
     ylim = range(c(els_run$gam, 1 / C), na.rm = TRUE))
abline(h = 1 / C, lty = 2, lwd = 2, col = "black")
legend("topright", legend = c(expression(gamma[n] ~ "(ELS)"), "1/C"),
       col = c("blue", "black"), lty = c(1, 2), lwd = 2, bty = "n")
dev.off()

# ---- Task 4 numbers for the discussion ----
# gamma_n = (1/2) / Rayleigh quotient of G at g, and lam_min(G) <= RQ <= lam_max(G) = C/2,
# so gamma_n >= 1/C always, with equality only when g is the top eigenvector.
cat("ELS step sizes: min =", min(els_run$gam, na.rm = TRUE),
    " max =", max(els_run$gam, na.rm = TRUE), " 1/C =", 1 / C, "\n")
cat("Ratio gamma_n * C (first 10):", round(els_run$gam[1:10] * C, 3), "\n")
cat("Upper bound 1/lam_min_pos =", 1 / lam_min_pos, "\n")
