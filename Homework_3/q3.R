
# -----------------------------------------------------------------------------
# PROBLEM 3 - Auto-MPG standardization and correlation (Dataset B)
# -----------------------------------------------------------------------------
# 1. Train = first 300 permuted rows, test = next 92.
# 2. Arrange as TRI (7 x 300), TRO (300 x 1), TEI (7 x 92), TEO (92 x 1).
#    Feature order: cylinders, displacement, horsepower, weight,
#                   acceleration, model_year, origin.
# 3. Standardize each feature row of TRI with training mean / population sd
#    (n = 300); apply the same transform to TEI. Do NOT touch TRO/TEO.
# Task 3.1: 8 x 8 absolute Pearson correlation matrix (7 features + mpg).
#           Which input correlates most with mpg?
# Task 3.2: Fig. 5 - TRO vs cylinders. Fig. 6 - TRO vs weight.
#           Fig. 7 - 3x3 grid of all 7 features vs mpg; jitter cylinders,
#           model_year, origin.

