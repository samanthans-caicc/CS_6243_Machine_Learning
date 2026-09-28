# Problem 5 (18 points): Low-Dimensional Geometry and Visualization

#   1. Start with if (!exists("sigma", inherits = FALSE)) source("q2.R") so you inherit U, X_centered, X, mu, labels,
#      original_indices, and r.

# inherits = FALSE matters: stats::sigma() is a real function, so a bare
# exists("sigma") is always TRUE and this guard would never fire.
if (!exists("sigma", inherits = FALSE)) source("q2.R")

#   2. Form Q2 <- U[, 1:2] and project with P2 <- crossprod(Q2, X_centered), giving a 2 x N matrix.
Q2 <- U[, 1:2]
P2 <- crossprod(Q2, X_centered)

#   3. For each class c in 0, 1, 9, compute the projected centroid as the row means of the columns of P2 whose
#      label is c, giving a length-2 vector per class.
classes <- c(0, 1, 9)
class_cols <- c("0" = "#D55E00", "1" = "#009E73", "9" = "#0072B2")

# rowMeans, not colMeans: P2 is 2 x N, so the samples are COLUMNS and the two
# coordinates are ROWS. colMeans here returns one value per image (N_c of them)
# instead of the length-2 centroid the PDF asks for.
centroids <- list()
for (c in classes) {
  centroids[[as.character(c)]] <- rowMeans(P2[, labels == c, drop = FALSE])
}

#   4. Report the three projected centroids as a small labeled table before plotting them.
table_centroids <- do.call(rbind, centroids)
colnames(table_centroids) <- c("PC1", "PC2")
cat("Projected class centroids:\n")
print(table_centroids, digits = 8)

png("fig5_projection.png", width = 1600, height = 1400, res = 200, type = "cairo")
plot(P2[1, ], P2[2, ],
     col = adjustcolor(class_cols[as.character(labels)], alpha.f = 0.25),
     pch = 19, cex = 0.35,
     xlab = "First principal direction",
     ylab = "Second principal direction",
     main = "MNIST digits 0, 1, 9 projected onto the first two left singular vectors")
for (c in classes) {
  # white halo underneath so the centroid marker stays visible over dense points
  points(centroids[[as.character(c)]][1], centroids[[as.character(c)]][2],
         col = "white", pch = 19, cex = 2.6)
  points(centroids[[as.character(c)]][1], centroids[[as.character(c)]][2],
         col = class_cols[as.character(c)], pch = 8, cex = 2.2, lwd = 3)
}
legend("topright",
       legend = c("digit 0", "digit 1", "digit 9", "class centroid"),
       col = c(class_cols, "black"), pch = c(19, 19, 19, 8), bty = "n")
dev.off()

#   5. Plot Fig. 5 as a scatter of P2[1, ] against P2[2, ], colored by class, and save it to its own PNG.

#   6. Use a small point size and semi-transparent colors, since 18,614 points will otherwise overplot into solid
#      blobs and hide the class structure.

#   7. Draw the three centroids on top with a distinct plotting symbol and a heavier size so they read clearly
#      against the scatter.

#   8. Add a legend naming the three digits and label both axes as the first and second principal directions.

#   9. Find the representative images by taking the first column of X whose label is 0, the first whose label is
#      1, and the first whose label is 9, using which(labels == c)[1].
# Name these rep_cols, NOT original_indices. q0.r already defines original_indices
# as the 0-based position of every retained image in the unfiltered 60,000-image
# split, and overwriting it would silently destroy the number the PDF asks you to
# report. rep_cols are 1-based COLUMN positions in the filtered matrix X.
rep_cols <- c(
    which(labels == 0)[1],
    which(labels == 1)[1],
    which(labels == 9)[1]
)
names(rep_cols) <- c("0", "1", "9")

original_images <- X[, rep_cols]

#  10. Report original_indices[j] for each of those three columns, which is the zero-based index in the unfiltered
#      60,000-image training split, and state plainly that this is not the post-filter index.

rep_table <- data.frame(
  digit            = c(0, 1, 9),
  mnist_index_0based = original_indices[rep_cols],  # report THIS one
  column_in_X        = rep_cols                     # not this one
)
cat("Representative images (report the zero-based MNIST index, not the column):\n")
print(rep_table, row.names = FALSE)
stopifnot(labels[rep_cols] == c(0, 1, 9))

#  10. Report original_indices[j] for each of those three columns, which is the zero-based index in the unfiltered
#      60,000-image training split, and state plainly that this is not the post-filter index.
    
#  11. Use K values 5, 20, and 100 as given. All three are below r = 632, so no substitution is needed, and you
#      should say so in the report since the PDF asks you to disclose any substitution.

#  12. Reconstruct each representative with x_hat <- mu + Q %*% crossprod(Q, x - mu), where Q <- U[, 1:K] and x is
#      the uncentered column X[, j].
reconstructions <- list()
for (K in c(5, 20, 100)) {
    Q <- U[, 1:K]
    x_hat <- sapply(1:ncol(original_images), function(j) {
        mu + Q %*% crossprod(Q, original_images[, j] - mu)
    })
    reconstructions[[as.character(K)]] <- x_hat
}

#  13. Compute all reconstructions before any display step, and do not clip the values, since the PDF forbids
#      clipping prior to numerical calculation.


#  14. Note that the 784-vectors unflatten transposed, because mnist.R stores each image transposed. Display with
#      image(matrix(v, 28, 28)[, 28:1], col = gray.colors(256)) so the digit appears upright.

#  15. Sanity-check that one original renders as a recognizable digit before trusting the whole figure, since a
#      wrong orientation is easy to miss in low-rank reconstructions.

#  16. Compute a single zlim across all twelve panels with range() over every original and reconstruction, then
#      pass that same zlim to every image() call so the grayscale mapping is common.
zlim <- range(c(original_images, unlist(reconstructions)))
cat("Common display limits across all 12 panels:", zlim, "\n")
cat("No clipping was applied at any stage, for display or otherwise.\n")

#  17. Plot Fig. 6 as a 3 x 4 grid using par(mfrow = c(3, 4)), one row per digit, columns being the original and
#      then K = 5, 20, 100.
# show_digit: one 784-vector as an upright 28 x 28 panel, on the shared zlim.
# The [, 28:1] flip is required because mnist.R stores each image transposed.
show_digit <- function(v, title) {
    image(matrix(v, 28, 28)[, 28:1], col = gray.colors(256), zlim = zlim,
          axes = FALSE, xlab = "", ylab = "", main = title)
    box()
}

png("fig6_reconstructions.png", width = 1600, height = 1300, res = 200, type = "cairo")
par(mfrow = c(3, 4), mar = c(1, 1, 2.5, 1))
for (i in seq_len(ncol(original_images))) {
    # Original FIRST, then increasing K, so the row reads left to right as the
    # rank grows. The caption has to match this order.
    show_digit(original_images[, i],
               paste0("digit ", names(rep_cols)[i], ": original"))
    for (K in c(5, 20, 100)) {
        show_digit(reconstructions[[as.character(K)]][, i], paste0("K = ", K))
    }
}
par(mfrow = c(1, 1))
dev.off()


#  18. If you clip anything for display only, say so explicitly in the caption, as the PDF requires that
#      disclosure.


#  19. In the geometry discussion, describe how much the three classes separate in two dimensions, which pair
#      overlaps most, and how that relates to E(2) from q2.R being small.

#  20. In the reconstruction discussion, describe what returns at each K: coarse blob and orientation at 5, digit
#      identity at 20, and stroke detail at 100.

#  21. For the proposed classifier, explain that it measures how far a test image sits from each class subspace
#      after removing that class mean, and assigns the class whose subspace explains the image best.

#  22. State that the residual is small when x lies near the class subspace, so the argmin picks the best-fitting
#      class, and note this is a per-class PCA model rather than a single global one.

#  23. Identify the training data it needs: labeled training images per class, from which you estimate mu_c and
#      Q_K,c by running an SVD on each class separately, plus a held-out set to choose K.

#  24. Say clearly that you are not required to implement or evaluate this classifier, so the report explains the
#      rationale only.
