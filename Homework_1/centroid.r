if (!exists("centroids")) source("q1.R")


# Show each centered centroid as a 28x28 image.
# Red = more ink than average, blue = less ink than average.
lim <- max(abs(unlist(centroids)))
pal <- colorRampPalette(c("blue", "white", "red"))(101)
png("centroids.png", width = 1500, height = 550, res = 150)
par(mfrow = c(1, 3), mar = c(1, 1, 2, 1))
for (c in classes) {
  v <- centroids[[as.character(c)]]
  image(matrix(v, 28, 28)[, 28:1], col = pal, zlim = c(-lim, lim),
        axes = FALSE, main = paste("Digit", c, "minus average"))
}
dev.off()

# A few summary numbers per class
for (c in classes) {
  v <- centroids[[as.character(c)]]
  cat("Class", c, ": min =", round(min(v), 2), " max =", round(max(v), 2),
      " pixels exactly 0 =", sum(v == 0), "\n")
}