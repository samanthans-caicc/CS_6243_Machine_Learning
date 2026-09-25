# Next step: write Homework_1/mnist.R with load_mnist(): open the file, skip 80 bytes, readBin the raw bytes, then shape the data. After
# that, run Rscript q0.r to see where it breaks next

load_mnist <- function() {
  # Open the file
  con <- file("../mnist/x_train.npy", "rb")
  # Skip the first 80 bytes (header)
  seek(con, 80)
  # Read the raw bytes and convert to integers 0..255
  raw_data <- as.integer(readBin(con, what = "raw", n = .Machine$integer.max))
  close(con)

  # Reshape into a 3D array: 28 x 28 x num_images.
  # R fills the first index fastest, but the file stores each image row by row,
  # so images[, , i] is the image TRANSPOSED. as.numeric(images[, , i]) is
  # therefore already the row-major vector that image.reshape(-1) gives in Python.
  images <- array(raw_data, dim = c(28, 28, length(raw_data) / (28 * 28)))

  # Load the labels (.npy header is also 80 bytes)
  con <- file("../mnist/y_train.npy", "rb")
  seek(con, 80)
  raw_labels <- readBin(con, what = "raw", n = .Machine$integer.max)
  close(con)
  labels <- as.integer(raw_labels)  # digits 0-9; do NOT subtract 1

  list(images = images, labels = labels)
}