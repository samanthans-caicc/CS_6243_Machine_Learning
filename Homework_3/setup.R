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
# Linux: use the venv that has numpy + scikit-learn (must come before import())
use_condaenv("base", required = TRUE)
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
# The file is WHITESPACE-delimited (not comma), so use read.table, not read.csv.
# car_name is in double quotes, which read.table handles by default.
auto <- read.table(
  "https://archive.ics.uci.edu/ml/machine-learning-databases/auto-mpg/auto-mpg.data",
  header = FALSE,
  na.strings = "?",
  col.names = c("mpg", "cylinders", "displacement", "horsepower", "weight",
                "acceleration", "model_year", "origin", "car_name")
)
auto <- auto[!is.na(auto$horsepower), ]
rownames(auto) <- NULL
stopifnot(nrow(auto) == 392)

# 5. Permutations with numpy so they match numpy.random.seed(42) exactly.
#    Seeds/sizes are passed as integers (42L) - numpy rejects R doubles here.
#    numpy returns 0-based indices, so add 1 for R.
#    We re-seed with 42 before each permutation (state this in the report).
np$random$seed(42L)
perm <- as.integer(np$random$permutation(length(x))) + 1L        # California, N = 20640

np$random$seed(42L)
perm_auto <- as.integer(np$random$permutation(nrow(auto))) + 1L  # Auto-MPG, N = 392

stopifnot(length(perm) == 20640, length(perm_auto) == 392)