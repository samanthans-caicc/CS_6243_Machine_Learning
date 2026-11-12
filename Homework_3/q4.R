# -----------------------------------------------------------------------------
# PROBLEM 4 - Univariate polynomial regression on each Auto-MPG feature
# -----------------------------------------------------------------------------
# 1. Double loop: 7 features x M in {0,1,2,3,4,5}, reusing fit_poly.
#    Store test MSEs in a 7 x 6 matrix.
#    - M = 0 is bias only -> make sure a 1-column matrix doesn't collapse
#      to a vector.
#    - Sanity check: all features give the same MSE at M = 0.
# Task 4.1: Fig. 8 - weight, M = 3. Dense raw grid over [1613, 5140].
#           Apply TRAINING x stats, then TRAINING column stats, predict,
#           plot curve vs RAW grid with TRO vs weight scatter.
#           Labels, legend, title.
# Task 4.2: Fig. 9 - same for horsepower, grid over [46, 230].
# Task 4.3: Fig. 10 - test MSE vs M, all 7 features, different colors.
# Task 4.4: Discussion - most predictive input; relate to Problem 3
#           correlations.

