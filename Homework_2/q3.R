# Problem 3 (25 points): GD with Fixed Step Size
# This file also defines the shared GD engine and timing helper that
# q4.R, q5.R and q6.R reuse.

if (!exists("L_opt")) source("q2.R")

# ---------------------------------------------------------------------------
# Shared GD engine
#   method = "fss" : gamma_n = gamma (fixed)
#   method = "els" : gamma_n = ||g||^2 / (2 g^T G g)   (closed-form exact LS step)
#   method = "bls" : start at gamma = 1, shrink by beta until Armijo holds
#
# Returns iterates X (K x 31, column n+1 holds x_n), accepted steps gam
# (NA when not applicable), updates performed, objective-evaluation count,
# backtracking reductions, and a status string.
#
# Objective-evaluation counting (BLS only; FSS and ELS need none):
#   * L(x_0) is evaluated once at the start (1 evaluation).
#   * Each Armijo trial L(x - gamma g) costs 1 evaluation, rejected or accepted.
#   * The accepted trial value is cached and reused as L(x_n) for the next
#     iteration, so it is never recomputed.
# ---------------------------------------------------------------------------
gd_run <- function(method, gamma = NULL, beta = NULL, c_arm = NULL,
                   max_backtracks = 100, curv_tol = 1e-300) {
  x   <- x_init
  X   <- matrix(NA_real_, K, n_iter + 1)
  X[, 1] <- x
  gam <- rep(NA_real_, n_iter)
  n_evals <- 0
  n_back  <- 0
  n_done  <- 0
  status  <- "completed 30 updates"

  if (method == "bls") {
    L_cur <- L_fun(x); n_evals <- n_evals + 1
  }

  for (n in 1:n_iter) {
    g  <- 2 * drop(G %*% x - v)
    gg <- sum(g^2)

    # Stopping check before computing a step
    if (sqrt(gg) <= grad_tol) {
      status <- sprintf("gradient tolerance met at x_%d", n - 1)
      break
    }

    if (method == "fss") {
      step <- gamma

    } else if (method == "els") {
      Gg   <- drop(G %*% g)
      curv <- sum(g * Gg)
      # Safeguard: g^T G g = ||Ag||^2 ~ 0 means no curvature along g, the
      # closed-form step is undefined. Stop and report instead of guessing.
      if (curv <= curv_tol) {
        status <- sprintf("ELS degenerate curvature at x_%d", n - 1)
        break
      }
      step <- gg / (2 * curv)

    } else if (method == "bls") {
      step <- 1
      accepted <- FALSE
      for (b in 0:max_backtracks) {
        L_trial <- L_fun(x - step * g); n_evals <- n_evals + 1
        if (L_trial <= L_cur - c_arm * step * gg) { accepted <- TRUE; break }
        step <- beta * step
        n_back <- n_back + 1
      }
      # Safeguard: if Armijo never holds, stop and report the failure rather
      # than silently taking some arbitrary step.
      if (!accepted) {
        status <- sprintf("BLS failed after %d reductions at x_%d", max_backtracks, n - 1)
        break
      }
      L_cur <- L_trial   # cache: becomes L(x_n)
    }

    x <- x - step * g
    X[, n + 1] <- x
    gam[n] <- step
    n_done <- n
  }

  # Early stop: copy x_m forward only to complete the plots; gamma stays NA
  if (n_done < n_iter) {
    for (j in (n_done + 2):(n_iter + 1)) X[, j] <- X[, n_done + 1]
  }

  list(X = X, gam = gam, n_updates = n_done, n_evals = n_evals,
       n_backtracks = n_back, status = status)
}

# Diagnostics computed AFTER the run from saved iterates (not timed, not counted)
obj_path  <- function(run) apply(run$X[, -1, drop = FALSE], 2, L_fun)       # L(x_1..x_30)
grad_path <- function(run) apply(run$X, 2, function(x) sqrt(sum(g_fun(x)^2))) # ||g(x_0..x_30)||

# Wall-clock timing: one untimed warm-up, then n_rep timed runs.
# Sys.time() has microsecond resolution on Linux (R's analogue of perf_counter).
n_rep <- 21
time_runs <- function(run_fun, n_rep_ = n_rep) {
  invisible(run_fun())                      # warm-up
  t <- numeric(n_rep_)
  for (i in 1:n_rep_) {
    t0 <- Sys.time()
    invisible(run_fun())
    t[i] <- as.numeric(Sys.time() - t0, units = "secs")
  }
  t
}

# ---------------------------------------------------------------------------
# Task 1: run FSS for each p and plot Fig. 1
# ---------------------------------------------------------------------------
p_set <- c(0.1, 0.5, 1, 1.5, 2)
fss_runs <- lapply(p_set, function(p) gd_run("fss", gamma = p / C))
names(fss_runs) <- paste0("p=", p_set)
fss_obj <- sapply(fss_runs, obj_path)       # 30 x 5

cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
png("fig1_fss.png", width = 1600, height = 1100, res = 200)
matplot(1:n_iter, fss_obj, type = "l", lty = 1, lwd = 2, col = cols, log = "y",
        xlab = "Iteration n", ylab = expression(L(x[n]) ~ "(log scale)"),
        main = "Fig. 1: GD with fixed step size gamma = p / C")
abline(h = L_opt, lty = 2, lwd = 2, col = "black")
legend("topright", legend = c(paste("p =", p_set), expression(L[opt])),
       col = c(cols, "black"), lty = c(rep(1, 5), 2), lwd = 2, bty = "n")
dev.off()

# ---------------------------------------------------------------------------
# Task 2: final gaps and iterate behaviour
# ---------------------------------------------------------------------------
fss_table <- data.frame(
  p           = p_set,
  gamma       = p_set / C,
  updates     = sapply(fss_runs, `[[`, "n_updates"),
  final_gap   = fss_obj[n_iter, ] - L_opt,
  dist_to_xopt = sapply(fss_runs, function(r) sqrt(sum((r$X[, n_iter + 1] - x_opt)^2))),
  grad_norm_30 = sapply(fss_runs, function(r) tail(grad_path(r), 1))
)
cat("\nProblem 3 summary:\n"); print(fss_table, digits = 6)

# Look at the iterates for p = 2: the component along the top eigenvector of H
# gets multiplied by (1 - 2) = -1 every step, so it flips sign but never shrinks.
eH <- eigen(H, symmetric = TRUE)
q1 <- eH$vectors[, 1]
proj_top <- function(run) drop(crossprod(q1, run$X - x_opt))   # n = 0..30
cat("\nComponent of x_n - x_opt along the top eigenvector of H (n = 0..6):\n")
for (nm in names(fss_runs)) cat(sprintf("%-6s", nm), round(proj_top(fss_runs[[nm]])[1:7], 6), "\n")

# Optional supporting plot of that oscillation
png("fig1b_fss_top_component.png", width = 1600, height = 1100, res = 200)
matplot(0:n_iter, sapply(fss_runs, proj_top), type = "o", pch = 20, lty = 1,
        col = cols, xlab = "Iteration n", ylab = "q1^T (x_n - x_opt)",
        main = "Error along top eigenvector of H")
abline(h = 0, lty = 3)
legend("topright", legend = paste("p =", p_set), col = cols, lty = 1, pch = 20, bty = "n")
dev.off()

# ---------------------------------------------------------------------------
# Task 3: theory. Error contracts in every positive-eigenvalue direction iff
# |1 - gamma * lambda_i| < 1 for all lambda_i > 0  <=>  0 < gamma < 2 / lambda_max = 2 / C,
# i.e. 0 < p < 2. The slowest factor is max(|1 - gamma lam_min_pos|, |1 - gamma C|).
# ---------------------------------------------------------------------------
rate <- sapply(p_set, function(p) max(abs(1 - (p / C) * lam_pos)))
cat("\nGuaranteed range: 0 < gamma < 2/C =", 2 / C, " (0 < p < 2)\n")
cat("Worst-case contraction factor max_i |1 - gamma lambda_i| per p:\n")
print(data.frame(p = p_set, factor = rate))
cat("Null-space dimension of H:", K - r,
    "(components there are multiplied by 1 every step and never change)\n")
