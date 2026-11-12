# -----------------------------------------------------------------------------
# SETUP (do once)
# -----------------------------------------------------------------------------
# 1. Load reticulate and import numpy (needed so the permutations match
#    numpy.random.seed(42) exactly - R's set.seed(42) gives different splits).
# 2. Record R version + package versions (sessionInfo(), py_config()) for the
#    Reproducibility section.
# 3. Load Dataset A (California Housing) via sklearn through reticulate.
#    x = column 1 (MedInc), y = target (units of $100k).
# 4. Load Dataset B (Auto-MPG) from the UCI URL.
#    - treat "?" as missing
#    - drop rows with missing horsepower
#    - sanity check: 392 rows remain
# 5. Seed numpy with 42, make the California permutation (N = 20640).
#    Re-seed with 42, make the Auto-MPG permutation (N = 392).
#    Add 1 to both (R is 1-indexed). Save them; never regenerate.
#    -> State the re-seeding choice in the report.

library(reticulate)
np <- import("numpy")

# 2. Record R version + package versions
sessionInfo()
py_config()

# 3. Load Dataset A (California Housing)
sklearn <- import("sklearn.datasets")
california <- sklearn$fetch_california_housing()
x <- california$data[, 1]
y <- california$target

# 4. Load Dataset B (Auto-MPG)
auto <- read.csv(
  "https://archive.ics.uci.edu/ml/machine-learning-databases/auto-mpg/auto-mpg.data",
  header = FALSE,
  na.strings = "?"
)
auto <- auto[!is.na(auto$V4), ]
stopifnot(nrow(auto) == 392)

