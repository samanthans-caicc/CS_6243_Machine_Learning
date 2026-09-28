# Problem 7 (18 points): SVD Analysis of Linear Systems

#   1. This problem is self-contained. Do not source q0.r through q6.R; none of the MNIST objects are used here.

#   2. Type in the four matrices A, B, C, D and the six right-hand sides yA, yB, yB2, yC, yC2, yD, yD2 exactly as
#      printed, keeping all eight decimals.

A <- matrix(c(
  -2.74125009,  2.24215689, -0.60553211, -0.16755625,
  -0.34868395,  0.29538923, -0.45259498,  0.50015934,
   2.49664208,  0.27798324,  2.00739274,  0.21978030
), nrow = 3, ncol = 4, byrow = TRUE)

yA <- c(0.61339829, 0.11012282, -0.06426754)

B <- matrix(c(
  -2.74125009, -0.34868395, 2.49664208,
   2.24215689,  0.29538923, 0.27798324,
  -0.60553211, -0.45259498, 2.00739274,
  -0.16755625,  0.50015934, 0.21978030
), nrow = 4, ncol = 3, byrow = TRUE)

yB  <- c(0.66761214,  0.35931116, 0.74289966, 0.02979187)
yB2 <- c(0.24982762, -0.45768269, 0.22778277, 0.63413920)

C <- matrix(c(
  0.31997336, 0.43316234, -0.33457014, -0.34017903,
  1.12969075, 1.52931319, -1.18122581, -1.20102843,
  0.20087760, 0.27193705, -0.21004138, -0.21356262
), nrow = 3, ncol = 4, byrow = TRUE)

yC  <- c( 0.14216640, 0.50192948,  0.08925132)
yC2 <- c(-1.01480112, 0.41152110, -0.45229071)

D <- matrix(c(
   0.07999334,  0.28242269,  0.05021940,
   0.10829058,  0.38232830,  0.06798426,
  -0.08364254, -0.29530645, -0.05251035,
  -0.08504476, -0.30025711, -0.05339065
), nrow = 4, ncol = 3, byrow = TRUE)

yD  <- c(0.41615372,  0.56336601, -0.43513813, -0.44243299)
yD2 <- c(0.47277025, -0.64357627,  1.30059591,  1.42694800)

#   3. Build each matrix with byrow = TRUE and then check dim() against the PDF, since a silent transpose here
#      corrupts every downstream answer.

#   4. Confirm A is 3 x 4, B is 4 x 3, C is 3 x 4, and D is 4 x 3 before going further.

stopifnot(dim(A) == c(3, 4), dim(B) == c(4, 3), dim(C) == c(3, 4), dim(D) == c(4, 3))
stopifnot(length(yA) == 3, length(yB) == 4, length(yB2) == 4)
stopifnot(length(yC) == 3, length(yC2) == 3, length(yD) == 4, length(yD2) == 4)

# Transcription check: the PDF prints B as the transpose of A, so this must hold exactly.
# It catches a mistyped digit in either matrix, which is the likeliest error in this problem.
stopifnot(identical(B, t(A)))

cat("dim(A):", dim(A), " dim(B):", dim(B), " dim(C):", dim(C), " dim(D):", dim(D), "\n")
cat("B == t(A):", identical(B, t(A)), "\n")

#   5. Set up the seven cases as a list of (G, f, name) triples so one function handles all of them and the report
#      stays consistent.

cases <- list(
  list(G = A, f = yA, name = "A"),
  list(G = B, f = yB, name = "B"),
  list(G = B, f = yB2, name = "B2"),
  list(G = C, f = yC, name = "C"),
  list(G = C, f = yC2, name = "C2"),
  list(G = D, f = yD, name = "D"),
  list(G = D, f = yD2, name = "D2")
)

#   6. Set options(digits = 10) or format every printed value explicitly, since the PDF demands at least six
#      decimal places.

options(digits = 10)

#   7. Write one function that takes G and f and returns everything the PDF asks for, then call it seven times.

analyze_case <- function(G, f, name) {
  s <- svd(G, nu = nrow(G), nv = ncol(G))
  tauG <- 1e-7 * s$d[1]
  rG <- sum(s$d > tauG)
  m <- nrow(G)
  n <- ncol(G)
  
  cat("Case:", name, "\n")
  cat("m x n:", m, "x", n, "\n")
  cat("Singular values:", s$d, "\n")
  cat("tauG:", tauG, "\n")
  cat("rG:", rG, "\n")
  
  shape <- if (m > n) {
    "overdetermined"
  } else if (m == n) {
    "square"
  } else {
    "underdetermined"
  }
  cat("System is", shape, "\n")

  # Step 25: flag rank deficiency, which is what drives the interesting cases.
  cat("min(m, n):", min(m, n),
      if (rG < min(m, n)) "-> rank deficient\n" else "-> full rank\n")

  Ur <- s$u[, 1:rG, drop = FALSE]
  P <- tcrossprod(Ur)
  # print(), not cat(): cat flattens P to one row of numbers, and report item 3
  # asks for the numerical entries as an actual matrix.
  cat("Projector P onto Range(G):\n")
  print(P)
  cat("Verification of projector P:\n")
  cat("norm(P - t(P), 'F'):", norm(P - t(P), "F"), "\n")
  cat("norm(P %*% P - P, 'F'):", norm(P %*% P - P, "F"), "\n")

  Pf <- as.vector(P %*% f)
  rf <- f - Pf
  rf_norm <- sqrt(sum(rf^2))
  cat("Pf:", Pf, "\n")
  cat("Residual rf:", rf, "\n")
  cat("norm(rf, '2'):", rf_norm, "\n")

  threshold <- 1e-7 * max(1, sqrt(sum(f^2)))
  in_range <- rf_norm <= threshold
  cat("threshold:", threshold, "\n")
  if (in_range) {
    cat("f is in Range(G)\n")
  } else {
    cat("f is not in Range(G)\n")
  }

  # Report item 5: consistency and rank together decide the solution count.
  verdict <- if (!in_range) {
    "no exact solution"
  } else if (rG == n) {
    "a unique exact solution"
  } else {
    "infinitely many exact solutions"
  }
  cat("Conclusion: the system has", verdict, "\n")
  cat("  justification:", if (in_range) "consistent" else "inconsistent",
      "with rG =", rG, "and n =", n, "\n")

  # drop = FALSE on s$v as well, so a rank-1 case stays a matrix instead of
  # falling back on R vector-matrix conformability rules.
  Gpinv <- s$v[, 1:rG, drop = FALSE] %*% diag(1/s$d[1:rG], rG, rG) %*% t(Ur)
  w <- as.vector(Gpinv %*% f)
  res_norm <- sqrt(sum((G %*% w - f)^2))

  # Report items 6 and 7: the same w, named for what it is in each case.
  w_label <- if (!in_range) {
    "minimum-norm least-squares solution wLS"
  } else if (rG == n) {
    "unique exact solution w"
  } else {
    "minimum-l2-norm solution w_dagger"
  }
  cat(w_label, ":\n")
  cat(w, "\n")
  cat("norm(w, '2'):", sqrt(sum(w^2)), "\n")
  cat("norm(G %*% w - f, '2'):", res_norm, "\n")
  if (!in_range) {
    # Item 7 cross-check: the least-squares residual must equal ||rf||.
    cat("  agrees with norm(rf, '2') to:", abs(res_norm - rf_norm), "\n")
  }

  invisible(list(
    case = name, m = m, n = n, shape = shape, tauG = tauG, rG = rG,
    rf_norm = rf_norm, in_range = in_range, verdict = verdict,
    w_norm = sqrt(sum(w^2)), res_norm = res_norm
  ))
}
#   8. Inside it, compute s <- svd(G, nu = nrow(G), nv = ncol(G)) so you have the full factors available.

results <- lapply(cases, function(cs) analyze_case(cs$G, cs$f, cs$name))

#  27. Present the seven cases as one table for the scalar quantities, with the projector matrices shown
#      separately above, so the report stays readable.

summary_table <- do.call(rbind, lapply(results, function(z) data.frame(
  case      = z$case,
  m         = z$m,
  n         = z$n,
  shape     = z$shape,
  rG        = z$rG,
  tauG      = z$tauG,
  resid     = z$rf_norm,
  in_range  = z$in_range,
  verdict   = z$verdict,
  norm_w    = z$w_norm,
  stringsAsFactors = FALSE
)))

cat("\n===================== SUMMARY OF ALL SEVEN CASES =====================\n")
print(summary_table, digits = 10, row.names = FALSE)

#   9. Compute tauG <- 1e-7 * s$d[1] and rG <- sum(s$d > tauG), using this problem-specific tolerance and not the
#      float64 tolerance from Problem 2.


#  10. Report m x n, all singular values of G, tauG, and rG for every case.
#  11. Classify the system as overdetermined when m > n, square when m = n, and underdetermined when m < n, and
#      report that label.

#  12. Form Ur <- s$u[, 1:rG, drop = FALSE] and P <- tcrossprod(Ur), the orthogonal projector onto Range(G).

#  13. Print the numerical entries of P as an actual matrix. The PDF forbids replacing it with a symbolic
#      expression.

#  14. Verify the projector by reporting norm(P - t(P), "F") and norm(P %*% P - P, "F"), both of which should land
#      near machine precision.

#  15. Compute Pf <- P %*% f and the residual rf <- f - Pf, then report both along with norm(rf, "2").

#  16. Use norm(as.matrix(rf), "F") or sqrt(sum(rf^2)) for the Euclidean norm, since norm() on a plain R vector
#      does not do what you expect.

#  17. Declare f in Range(G) if and only if the residual norm is at most 1e-7 * max(1, sqrt(sum(f^2))), using the
#      exact criterion the PDF states.

#  18. Build the truncated pseudoinverse yourself as Gpinv <- s$v[, 1:rG] %*% diag(1/s$d[1:rG], rG, rG) %*%
#      t(Ur), rather than calling a library pseudoinverse.

#  19. Pass rG twice to diag() so a rank-1 case yields a 1 x 1 matrix instead of collapsing to a scalar.

#  20. State in the report that MASS::ginv was not used, or that its tolerance was set to match, because the PDF
#      forbids a default library cutoff.

#  21. Never form t(G) %*% G anywhere. The PDF rules out normal equations as the primary method.

#  22. Decide the solution count from consistency and rank together: inconsistent means no exact solution,
#      consistent with rG = n means a unique one, and consistent with rG < n means infinitely many.

#  23. Report w <- Gpinv %*% f in every case, and label it the minimum-l2-norm solution when there are infinitely
#      many and the minimum-norm least-squares solution when there are none.

#  24. Report norm(G %*% w - f, "2") whenever there is no exact solution, and confirm it matches norm(rf, "2").

#  25. Expect the two rank-deficient matrices to drive the interesting cases, so check rG against min(m, n)
#      before writing the justification for each case.

#  26. Note that the paired right-hand sides are the point of the problem: one of each pair is built to lie in
#      Range(G) and the other is not, so expect different verdicts within the same matrix.

#  27. Present the seven cases as one table for the scalar quantities, with the projector matrices shown
#      separately, so the report stays readable.
